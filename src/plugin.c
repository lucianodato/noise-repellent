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
* \file plugin.c
* \author Luciano Dato
* \brief The main file for host interaction
*/

#include "plugin_state.c"
#include "stft_processor.c"

///---------------------------------------------------------------------

/**
* Enumeration of LV2 ports.
*/
typedef enum {
	NOISEREPELLENT_AMOUNT = 0,
	NOISEREPELLENT_NOFFSET = 1,
	NOISEREPELLENT_RELEASE = 2,
	NOISEREPELLENT_MASKING = 3,
	NOISEREPELLENT_T_PROTECT = 4,
	NOISEREPELLENT_WHITENING = 5,
	NOISEREPELLENT_N_LEARN = 6,
	NOISEREPELLENT_N_ADAPTIVE = 7,
	NOISEREPELLENT_RESET = 8,
	NOISEREPELLENT_RESIDUAL_LISTEN = 9,
	NOISEREPELLENT_ENABLE = 10,
	NOISEREPELLENT_LATENCY = 11,
	NOISEREPELLENT_INPUT = 12,
	NOISEREPELLENT_OUTPUT = 13,
} PIndex;

/**
* Struct for noise repellent instance, the host is going to use.
*/
typedef struct
{
	const float *input; //input of samples from host (changing size)
	float *output;		//output of samples to host (changing size)
	float sample_rate;  //Sample rate received from the host

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

	//STFT processing instance
	STFTprocessor *stft_processor; //The stft transform object

	//Plugin state instance
	Pstate *plugin_state;
} Nrepellent;

/**
* Instantiates the plugin.
*/
static LV2_Handle
instantiate(const LV2_Descriptor *descriptor, double rate, const char *bundle_path,
			const LV2_Feature *const *features)
{
	//Actual struct declaration
	Nrepellent *self = (Nrepellent *)calloc(1, sizeof(Nrepellent));

	// //Plugin state initialization
	// if (!ps_init(self->plugin_state, features))
	// {
	// 	//bail out: host does not support urid:map
	//     free(self);
	//     return NULL;
	// }

	//Sampling related
	self->sample_rate = (float)rate;

	//STFT related
	self->stft_processor = stft_p_init(self->sample_rate);

	return (LV2_Handle)self;
}

/**
* Used by the host to connect the ports of this plugin.
*/
static void
connect_port(LV2_Handle instance, uint32_t port, void *data)
{
	Nrepellent *self = (Nrepellent *)instance;

	switch ((PIndex)port)
	{
	case NOISEREPELLENT_AMOUNT:
		self->amount_of_reduction = (float *)data;
		break;
	case NOISEREPELLENT_NOFFSET:
		self->noise_thresholds_offset = (float *)data;
		break;
	case NOISEREPELLENT_RELEASE:
		self->release = (float *)data;
		break;
	case NOISEREPELLENT_MASKING:
		self->masking = (float *)data;
		break;
	case NOISEREPELLENT_WHITENING:
		self->whitening_factor_pc = (float *)data;
		break;
	case NOISEREPELLENT_N_LEARN:
		self->noise_learn_state = (float *)data;
		break;
	case NOISEREPELLENT_N_ADAPTIVE:
		self->adaptive_state = (float *)data;
		break;
	case NOISEREPELLENT_RESIDUAL_LISTEN:
		self->residual_listen = (float *)data;
		break;
	case NOISEREPELLENT_T_PROTECT:
		self->transient_protection = (float *)data;
		break;
	case NOISEREPELLENT_RESET:
		self->reset_profile = (float *)data;
		break;
	case NOISEREPELLENT_ENABLE:
		self->enable = (float *)data;
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

/**
* Main process function of the plugin.
*/
static void
run(LV2_Handle instance, uint32_t n_samples)
{
	Nrepellent *self = (Nrepellent *)instance;

	//Inform latency at run call
	*(self->report_latency) = (float)stft_p_get_latency(self->stft_processor);

	//Run the stft transform and process samples
	stft_p_run(self->stft_processor, n_samples, self->input, self->output, (int)*self->enable);
}

/**
* Cleanup and freeing memory.
*/
static void
cleanup(LV2_Handle instance)
{
	Nrepellent *self = (Nrepellent *)instance;

	stft_p_free(self->stft_processor);
	free(instance);
}

/**
* State saving of the noise profile.
*/
static LV2_State_Status
savestate(LV2_Handle instance, LV2_State_Store_Function store, LV2_State_Handle handle,
		  uint32_t flags, const LV2_Feature *const *features)
{
	//Nrepellent *self = (Nrepellent *)instance;

	// ps_savestate(self->plugin_state, store, handle, self->stft_processor->fft_size,
	// 			 self->stft_processor->fft_processor->fft_denoiser->noise_estimation->noise_window_count,
	// 			 self->stft_processor->fft_processor->fft_denoiser->noise_estimation->noise_profile);

	return LV2_STATE_SUCCESS;
}

/**
* State restoration of the noise profile.
*/
static LV2_State_Status
restorestate(LV2_Handle instance, LV2_State_Retrieve_Function retrieve,
			 LV2_State_Handle handle, uint32_t flags,
			 const LV2_Feature *const *features)
{
	//Nrepellent *self = (Nrepellent *)instance;

	// if(!ps_restorestate(self->plugin_state, retrieve, handle,
	// 				self->stft_processor->fft_processor->fft_denoiser->noise_estimation->noise_profile,
	// 				self->stft_processor->fft_processor->fft_denoiser->noise_estimation->noise_window_count,
	// 				self->stft_processor->fft_size, self->stft_processor->fft_size_2))
	// {
	// 	return LV2_STATE_ERR_NO_PROPERTY;
	// }

	return LV2_STATE_SUCCESS;
}

/**
* extension for additional interfaces.
*/
static const void *
extension_data(const char *uri)
{
	static const LV2_State_Interface state = {savestate, restorestate};
	if (!strcmp(uri, LV2_STATE__interface))
	{
		return &state;
	}
	return NULL;
}

/**
* Descriptor for linking methods.
*/
static const LV2_Descriptor descriptor =
{
	NOISEREPELLENT_URI,
	instantiate,
	connect_port,
	NULL,
	run,
	NULL,
	cleanup,
	extension_data
};

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
