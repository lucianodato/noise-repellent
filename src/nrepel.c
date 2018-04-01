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

/**
* \file nrepel.c
* \author Luciano Dato
* \brief The main file for host interaction
*/

#include <fftw3.h>
#include "lv2/lv2plug.in/ns/lv2core/lv2.h"
#include "lv2/lv2plug.in/ns/ext/urid/urid.h"
#include "lv2/lv2plug.in/ns/ext/atom/atom.h"
#include "lv2/lv2plug.in/ns/ext/state/state.h"

#include "stft.c"

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"

//STFT default values
#define FFT_SIZE 2048	//Size of the fft transform
#define BLOCK_SIZE 2048	//Size of the block of samples
#define INPUT_WINDOW 3   //0 HANN 1 HAMMING 2 BLACKMAN 3 VORBIS Input windows for STFT algorithm
#define OUTPUT_WINDOW 3  //0 HANN 1 HAMMING 2 BLACKMAN 3 VORBIS Output windows for STFT algorithm
#define OVERLAP_FACTOR 4 //4 is 75% overlap Values bigger than 4 will rescale correctly (if Vorbis windows is not used)

///---------------------------------------------------------------------

/**
* Enumeration of LV2 ports.
*/
typedef enum {
	NREPEL_AMOUNT = 0,
	NREPEL_NOFFSET = 1,
	NREPEL_RELEASE = 2,
	NREPEL_MASKING = 3,
	NREPEL_T_PROTECT = 4,
	NREPEL_WHITENING = 5,
	NREPEL_N_LEARN = 6,
	NREPEL_N_ADAPTIVE = 7,
	NREPEL_RESET = 8,
	NREPEL_RESIDUAL_LISTEN = 9,
	NREPEL_ENABLE = 10,
	NREPEL_LATENCY = 11,
	NREPEL_INPUT = 12,
	NREPEL_OUTPUT = 13,
} PortIndex;

/**
* Struct for THE noise repellent instance, the host is going to use.
*/
typedef struct
{
	const float *input; //input of samples from host (changing size)
	float *output;		//output of samples to host (changing size)
	float samp_rate;	//Sample rate received from the host

	//Parameters for the algorithm (user input)
	float *amount_of_reduction;		//Amount of noise to reduce in dB
	float *noise_thresholds_offset; //This is to scale the noise profile (over subtraction factor)
	float *release;					//Release time
	float *masking;					//Masking scaling
	float *whitening_factor_pc;		//Whitening amount of the reduction percentage
	float *noise_learn_state;		//Learn Noise state (Manual-Off-Auto)
	float *adaptive_state;			//Autocapture switch
	float *reset_profile;			//Reset Noise switch
	float *residual_listen;			//For noise only listening
	float *transient_protection;	//Multiplier for thresholding onsets with rolling mean
	float *enable;					//For soft bypass (click free bypass)
	float *report_latency;			//Latency necessary

	//Parameters values and arrays for the STFT
	STFTtransform* transform;  //The stft transform object


	// //LV2 state URID (Save and restore noise profile)
	// LV2_URID_Map *map;
	// LV2_URID atom_Vector;
	// LV2_URID atom_Int;
	// LV2_URID atom_Float;
	// LV2_URID prop_fftsize;
	// LV2_URID prop_nwindow;
	// LV2_URID prop_FFTp2;
} Nrepel;

/**
* Instantiates the plugin.
*/
static LV2_Handle
instantiate(const LV2_Descriptor *descriptor, double rate, const char *bundle_path,
			const LV2_Feature *const *features)
{
	//Actual struct declaration
	Nrepel *self = (Nrepel *)calloc(1, sizeof(Nrepel));

	// //Retrieve the URID map callback, and needed URIDs
	// for (int i = 0; features[i]; ++i)
	// {
	// 	if (!strcmp(features[i]->URI, LV2_URID__map))
	// 	{
	// 		self->map = (LV2_URID_Map *)features[i]->data;
	// 	}
	// }
	// if (!self->map)
	// {
	// 	//bail out: host does not support urid:map
	// 	free(self);
	// 	return NULL;
	// }

	// //For lv2 state (noise profile saving)
	// self->atom_Vector = self->map->map(self->map->handle, LV2_ATOM__Vector);
	// self->atom_Int = self->map->map(self->map->handle, LV2_ATOM__Int);
	// self->atom_Float = self->map->map(self->map->handle, LV2_ATOM__Float);
	// self->prop_fftsize = self->map->map(self->map->handle, NREPEL_URI "#fftsize");
	// self->prop_nwindow = self->map->map(self->map->handle, NREPEL_URI "#nwindow");
	// self->prop_FFTp2 = self->map->map(self->map->handle, NREPEL_URI "#FFTp2");

	//Sampling related
	self->samp_rate = (float)rate;

	//STFT related
	self->transform = stft_init(BLOCK_SIZE, FFT_SIZE, INPUT_WINDOW, OUTPUT_WINDOW, OVERLAP_FACTOR);

	return (LV2_Handle)self;
}

/**
* Used by the host to connect the ports of this plugin.
*/
static void
connect_port(LV2_Handle instance, uint32_t port, void *data)
{
	Nrepel *self = (Nrepel *)instance;

	switch ((PortIndex)port)
	{
	case NREPEL_AMOUNT:
		self->amount_of_reduction = (float *)data;
		break;
	case NREPEL_NOFFSET:
		self->noise_thresholds_offset = (float *)data;
		break;
	case NREPEL_RELEASE:
		self->release = (float *)data;
		break;
	case NREPEL_MASKING:
		self->masking = (float *)data;
		break;
	case NREPEL_WHITENING:
		self->whitening_factor_pc = (float *)data;
		break;
	case NREPEL_N_LEARN:
		self->noise_learn_state = (float *)data;
		break;
	case NREPEL_N_ADAPTIVE:
		self->adaptive_state = (float *)data;
		break;
	case NREPEL_RESIDUAL_LISTEN:
		self->residual_listen = (float *)data;
		break;
	case NREPEL_T_PROTECT:
		self->transient_protection = (float *)data;
		break;
	case NREPEL_RESET:
		self->reset_profile = (float *)data;
		break;
	case NREPEL_ENABLE:
		self->enable = (float *)data;
		break;
	case NREPEL_LATENCY:
		self->report_latency = (float *)data;
		break;
	case NREPEL_INPUT:
		self->input = (const float *)data;
		break;
	case NREPEL_OUTPUT:
		self->output = (float *)data;
		break;
	}
}

/**
* Main process function of the plugin.
*/
static void
run(LV2_Handle instance, uint32_t n_samples)
{
	Nrepel *self = (Nrepel *)instance;

	//Inform latency at run call
	*(self->report_latency) = (float)get_latency(self->transform);

	//Run the stft transform and process samples
	stft_run(self->transform, n_samples, self->input, self->output);
}

/**
* Cleanup and freeing memory.
*/
static void
cleanup(LV2_Handle instance)
{
	stft_free(instance->transform);
	free(instance);
}

// /**
// * Noise Profile state.
// */
// typedef struct
// {
// 	uint32_t child_size;
// 	uint32_t child_type;
// 	float array[(FFT_SIZE/2) + 1];
// } NProfile;

// /**
// * State saving of the noise profile.
// */
// static LV2_State_Status
// savestate(LV2_Handle instance, LV2_State_Store_Function store, LV2_State_Handle handle,
// 		  uint32_t flags, const LV2_Feature *const *features)
// {
// 	Nrepel *self = (Nrepel *)instance;

// 	NProfile *vector = (NProfile *)malloc(sizeof(NProfile));

// 	vector->child_type = self->atom_Float;
// 	vector->child_size = sizeof(float);

// 	store(handle, self->prop_fftsize, &self->transform->fft_size, sizeof(int), self->atom_Int,
// 		  LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

// 	store(handle, self->prop_nwindow, &self->noise_window_count, sizeof(float),
// 		  self->atom_Float, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

// 	memcpy(vector->array, self->noise_thresholds_p2, sizeof(vector->array));

// 	store(handle, self->prop_FFTp2, (void *)vector, sizeof(NProfile),
// 		  self->atom_Vector, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

// 	return LV2_STATE_SUCCESS;
// }

// /**
// * State restoration of the noise profile.
// */
// static LV2_State_Status
// restorestate(LV2_Handle instance, LV2_State_Retrieve_Function retrieve,
// 			 LV2_State_Handle handle, uint32_t flags,
// 			 const LV2_Feature *const *features)
// {
// 	Nrepel *self = (Nrepel *)instance;
// 	size_t size;
// 	uint32_t type;
// 	uint32_t valflags;

// 	const int32_t *fftsize = retrieve(handle, self->prop_fftsize, &size, &type, &valflags);
// 	if (!fftsize || type != self->atom_Int || *fftsize != self->transform->fft_size_2)
// 	{
// 		return LV2_STATE_ERR_NO_PROPERTY;
// 	}

// 	const void *vecFFTp2 = retrieve(handle, self->prop_FFTp2, &size, &type, &valflags);
// 	if (!vecFFTp2 || size != sizeof(NProfile) || type != self->atom_Vector)
// 	{
// 		return LV2_STATE_ERR_NO_PROPERTY;
// 	}

// 	//Deactivate any denoising before loading any noise profile
// 	self->noise_thresholds_availables = false;

// 	//Copy to local variables
// 	memcpy(self->noise_thresholds_p2, (float *)LV2_ATOM_BODY(vecFFTp2), (self->transform->fft_size_2 + 1) * sizeof(float));

// 	const float *wincount = retrieve(handle, self->prop_nwindow, &size, &type, &valflags);
// 	if (fftsize && type == self->atom_Float)
// 	{
// 		self->noise_window_count = *wincount;
// 	}

// 	//Reactivate denoising with restored profile
// 	self->noise_thresholds_availables = true;

// 	return LV2_STATE_SUCCESS;
// }

// /**
// * extension for additional interfaces.
// */
// static const void *
// extension_data(const char *uri)
// {
// 	static const LV2_State_Interface state = {savestate, restorestate};
// 	if (!strcmp(uri, LV2_STATE__interface))
// 	{
// 		return &state;
// 	}
// 	return NULL;
// }

/**
* Descriptor for linking methods.
*/
static const LV2_Descriptor descriptor =
	{
		NREPEL_URI,
		instantiate,
		connect_port,
		NULL,
		run,
		NULL,
		cleanup/*,
		extension_data*/};

/**
* Symbol export using the descriptor above.
*/
LV2_SYMBOL_EXPORT
const LV2_Descriptor *
lv2_descriptor(uint32_t index)
{
	switch (index)
	{
	case 0:
		return &descriptor;
	default:
		return NULL;
	}
}
