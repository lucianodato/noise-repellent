/*
noise-repellent -- Noise Reduction LV2

Copyright 2016 Luciano Dato <lucianodato@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/
*/

#include "denoise_parameters.h"
#include "fft_denoiser.h"
#include "lv2/atom/atom.h"
#include "lv2/core/lv2.h"
#include "lv2/core/lv2_util.h"
#include "lv2/log/logger.h"
#include "lv2/state/state.h"
#include "lv2/urid/urid.h"
#include "noise_profile.h"
#include "stft_processor.h"

#include <stdlib.h>
#include <string.h>

#define NOISEREPELLENT_URI "https://github.com/lucianodato/noise-repellent"
#define FFT_SIZE 2048
#define OVERLAP_FACTOR 4

typedef struct
{
	LV2_URID atom_Int;
	LV2_URID atom_Vector;
	LV2_URID plugin;
	LV2_URID atom_URID;
} URIs;

typedef struct
{
	LV2_URID property_fft_size;
	LV2_URID property_saved_noise_profile;
} State;

static inline void
map_uris(LV2_URID_Map *map, URIs *uris)
{
	uris->plugin = map->map(map->handle, NOISEREPELLENT_URI);
	uris->atom_Int = map->map(map->handle, LV2_ATOM__Int);
	uris->atom_Vector = map->map(map->handle, LV2_ATOM__Vector);
	uris->atom_URID = map->map(map->handle, LV2_ATOM__URID);
}

static inline void
map_state(LV2_URID_Map *map, State *state)
{
	state->property_fft_size = map->map(map->handle, NOISEREPELLENT_URI "#fftsize");
	state->property_saved_noise_profile = map->map(map->handle, NOISEREPELLENT_URI "#savednoiseprofile");
}

typedef enum
{
	NOISEREPELLENT_AMOUNT = 0,
	NOISEREPELLENT_NOISE_OFFSET = 1,
	NOISEREPELLENT_RELEASE = 2,
	NOISEREPELLENT_MASKING = 3,
	NOISEREPELLENT_TRANSIENT_PROTECT = 4,
	NOISEREPELLENT_WHITENING = 5,
	NOISEREPELLENT_NOISE_LEARN = 6,
	NOISEREPELLENT_RESIDUAL_LISTEN = 7,
	NOISEREPELLENT_ENABLE = 8,
	NOISEREPELLENT_LATENCY = 9,
	NOISEREPELLENT_INPUT = 10,
	NOISEREPELLENT_OUTPUT = 11,
} PortIndex;

typedef struct
{
	const float *input;
	float *output;
	float sample_rate;
	float *report_latency;

	// Features
	LV2_URID_Map *map;
	LV2_Log_Logger log;
	URIs uris;
	State state;

	DenoiseParameters *denoise_parameters;
	NoiseProfile *noise_profile;
	FFTDenoiser *fft_denoiser;
	STFTProcessor *stft_processor;
} NoiseRepellent;

static LV2_Handle instantiate(const LV2_Descriptor *descriptor, double rate, const char *bundle_path,
							  const LV2_Feature *const *features)
{
	NoiseRepellent *self = (NoiseRepellent *)calloc(1, sizeof(NoiseRepellent));

	const char *missing = lv2_features_query(
		features,
		LV2_LOG__log, &self->log.log, false,
		LV2_URID__map, &self->map, true,
		NULL);

	lv2_log_logger_set_map(&self->log, self->map);

	if (missing)
	{
		lv2_log_error(&self->log, "Missing feature <%s>\n", missing);
		free(self);
		return NULL;
	}

	map_uris(self->map, &self->uris);
	map_state(self->map, &self->state);

	self->sample_rate = (float)rate;

	self->noise_profile = noise_profile_initialize(FFT_SIZE / 2 + 1);
	self->denoise_parameters = plugin_parameters_initialize();
	self->fft_denoiser = fft_denoiser_initialize(self->sample_rate, FFT_SIZE, OVERLAP_FACTOR);
	self->stft_processor = stft_processor_initialize(self->fft_denoiser, FFT_SIZE, OVERLAP_FACTOR);

	return (LV2_Handle)self;
}

static void connect_port(LV2_Handle instance, uint32_t port, void *data)
{
	NoiseRepellent *self = (NoiseRepellent *)instance;

	switch ((PortIndex)port)
	{
	case NOISEREPELLENT_AMOUNT:
		set_plugin_parameters(self->denoise_parameters, (float *)data, REDUCTION_AMOUNT);
		break;
	case NOISEREPELLENT_NOISE_OFFSET:
		set_plugin_parameters(self->denoise_parameters, (float *)data, NOISE_RESCALE);
		break;
	case NOISEREPELLENT_RELEASE:
		set_plugin_parameters(self->denoise_parameters, (float *)data, RELEASE);
		break;
	case NOISEREPELLENT_MASKING:
		set_plugin_parameters(self->denoise_parameters, (float *)data, MASKING);
		break;
	case NOISEREPELLENT_WHITENING:
		set_plugin_parameters(self->denoise_parameters, (float *)data, WHITENING_FACTOR);
		break;
	case NOISEREPELLENT_NOISE_LEARN:
		set_plugin_parameters(self->denoise_parameters, (float *)data, LEARN_NOISE);
		break;
	case NOISEREPELLENT_RESIDUAL_LISTEN:
		set_plugin_parameters(self->denoise_parameters, (float *)data, RESIDUAL_LISTEN);
		break;
	case NOISEREPELLENT_TRANSIENT_PROTECT:
		set_plugin_parameters(self->denoise_parameters, (float *)data, TRANSIENT_PROTECTION);
		break;
	case NOISEREPELLENT_ENABLE:
		set_plugin_parameters(self->denoise_parameters, (float *)data, ENABLE);
		break;
	case NOISEREPELLENT_LATENCY:
		self->report_latency = (float *)data;
		break;
	case NOISEREPELLENT_INPUT:
		self->input = (const float *)data;
		break;
	case NOISEREPELLENT_OUTPUT:
		self->output = (float *)data;
		break;
	}
}

static void run(LV2_Handle instance, uint32_t n_samples)
{
	NoiseRepellent *self = (NoiseRepellent *)instance;

	load_denoise_parameters(self->fft_denoiser, self->denoise_parameters);

	*(self->report_latency) = (float)stft_processor_get_latency(self->stft_processor);
	stft_processor_run(self->stft_processor, self->noise_profile, n_samples, self->input, self->output);
}

static void cleanup(LV2_Handle instance)
{
	NoiseRepellent *self = (NoiseRepellent *)instance;

	noise_profile_free(self->noise_profile);
	plugin_parameters_free(self->denoise_parameters);
	fft_denoiser_free(self->fft_denoiser);
	stft_processor_free(self->stft_processor);
	free(instance);
}

static LV2_State_Status savestate(LV2_Handle instance, LV2_State_Store_Function store, LV2_State_Handle handle,
								  uint32_t flags, const LV2_Feature *const *features)
{
	NoiseRepellent *self = (NoiseRepellent *)instance;

	const int fft_size = FFT_SIZE;
	const float *noise_profile = get_noise_profile(self->noise_profile);

	store(handle, self->state.property_fft_size, &fft_size, sizeof(int), self->uris.atom_Int,
		  LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	store(handle, self->state.property_saved_noise_profile, (void *)noise_profile, sizeof(float) * fft_size / 2 + 1,
		  self->uris.atom_Vector, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	return LV2_STATE_SUCCESS;
}

static LV2_State_Status restorestate(LV2_Handle instance, LV2_State_Retrieve_Function retrieve,
									 LV2_State_Handle handle, uint32_t flags,
									 const LV2_Feature *const *features)
{
	NoiseRepellent *self = (NoiseRepellent *)instance;

	size_t size;
	uint32_t type;
	uint32_t valflags;

	const int *fftsize = (const int *)retrieve(handle, self->state.property_fft_size, &size, &type, &valflags);
	if (!fftsize || type != self->uris.atom_Int)
	{
		return LV2_STATE_ERR_NO_PROPERTY;
	}

	if (!*fftsize)
		set_spectral_size(self->stft_processor, FFT_SIZE);
	else
		set_spectral_size(self->stft_processor, *fftsize);

	const void *saved_noise_profile = retrieve(handle, self->state.property_saved_noise_profile, &size, &type, &valflags);
	if (!saved_noise_profile || size != sizeof(float) * *fftsize / 2 + 1 || type != self->uris.atom_Vector)
	{
		return LV2_STATE_ERR_NO_PROPERTY;
	}

	set_noise_profile(self->noise_profile, (float *)LV2_ATOM_BODY(saved_noise_profile));

	return LV2_STATE_SUCCESS;
}

static const void *extension_data(const char *uri)
{
	static const LV2_State_Interface state = {savestate, restorestate};
	if (!strcmp(uri, LV2_STATE__interface))
	{
		return &state;
	}
	return NULL;
}

static const LV2_Descriptor descriptor =
	{
		NOISEREPELLENT_URI,
		instantiate,
		connect_port,
		NULL,
		run,
		NULL,
		cleanup,
		extension_data};

LV2_SYMBOL_EXPORT const LV2_Descriptor *lv2_descriptor(uint32_t index)
{
	switch (index)
	{
	case 0:
		return &descriptor;
	default:
		return NULL;
	}
}
