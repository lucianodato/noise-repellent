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

#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <time.h>
#include <stdio.h>

#include "spectral_processing.c"
#include "lv2/lv2plug.in/ns/lv2core/lv2.h"
#include "lv2/lv2plug.in/ns/ext/urid/urid.h"
#include "lv2/lv2plug.in/ns/ext/atom/atom.h"
#include "lv2/lv2plug.in/ns/ext/state/state.h"

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"

//STFT default values
#define FFT_SIZE 2048 //Size of the fft transform
#define INPUT_WINDOW 2 //0 HANN 1 HAMMING 2 BLACKMAN Input windows for STFT algorithm
#define OUTPUT_WINDOW 0 //0 HANN 1 HAMMING 2 BLACKMAN Input windows for STFT algorithm
#define FRAME_SIZE 2048 //Size of the analysis windows
#define OVERLAP_FACTOR 4 //4 is 75% overlap Values bigger than 4 will rescale correctly (2 won't work as is using windowing for analysis and Synthesis)
#define NOISE_ARRAY_STATE_MAX_SIZE 8192 //max alloc size of the noise_thresholds to save with the session. This will consider upto fs of 192 kHz fs


///---------------------------------------------------------------------

//LV2 CODE

typedef enum
{
	NREPEL_AMOUNT = 0,
	NREPEL_NOFFSET = 1,
	NREPEL_RELEASE = 2,
	NREPEL_MASKING = 3,
	NREPEL_A_CONTROL = 4,
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

//spectum struct for noise profile saving
typedef struct
{
	uint32_t child_size;
	uint32_t child_type;
	float array[NOISE_ARRAY_STATE_MAX_SIZE];
} FFTVector;

typedef struct
{
	const float* input; //input of samples from host (changing size)
	float* output; //output of samples to host (changing size)
	float samp_rate; //Sample rate received from the host

	//Parameters for the algorithm (user input)
	float* amount_of_reduction; //Amount of noise to reduce in dB
	float* noise_thresholds_offset; //This is to scale the noise profile (over subtraction factor)
	float* release; //Release time
	float* masking; //Masking scaling
	float* artifact_control; //Mixing between spectral gates and spectral subtraction
	float* whitening_factor_pc;	//Whitening amount of the reduction percentage
	float* noise_learn_state; //Learn Noise state (Manual-Off-Auto)
	float* adaptive_state; //Autocapture switch
	float* reset_profile; //Reset Noise switch
	float* residual_listen; //For noise only listening
	float* enable; //For soft bypass (click free bypass)
	float* report_latency; //Latency necessary

	//Control variables
	bool noise_thresholds_availables; //indicate whether a noise profile is available or no
	bool peak_detection;

	//Parameters values and arrays for the STFT
	int fft_size; //FFTW input size
	int fft_size_2; //FFTW half input size
	int window_option_input; //Type of input Window for the STFT
	int window_option_output; //Type of output Window for the STFT
	float overlap_factor; //oversampling factor for overlap calculations
	float overlap_scale_factor; //Scaling factor for conserving the final amplitude
	int hop; //Hop size for the STFT
	float* input_window; //Input Window values
	float* output_window; //Output Window values
	int frame_size; //Frame size
	int frame_fft_diff;	//Difference between the fft size and the window size
	float window_count; //Count windows for mean computing

	//Algorithm exta variables
	float tau; //time constant for soft bypass
	float wet_dry_target; //softbypass target for softbypass
	float wet_dry; //softbypass coeff
	float reduction_coeff; //Gain to apply to the residual noise
	float release_coeff; //Release coefficient for Envelopes
	float amount_of_reduction_linear;	//Reduction amount linear value
	float thresholds_offset_linear;	//Threshold offset linear value
	float whitening_factor; //Whitening amount of the reduction
	float prev_beta; //For the adaptive smoothing

	//Buffers for processing and outputting
	int input_latency;
	float* in_fifo; //internal input buffer
	float* out_fifo; //internal output buffer
	float* output_accum; //FFT output accumulator
	int read_ptr; //buffers read pointer

	//FFTW related arrays
	float* input_fft_buffer;
	float* output_fft_buffer;
	fftwf_plan forward;
	fftwf_plan backward;

	//Postfilter related
	float* input_fft_buffer_ps;
	float* output_fft_buffer_ps;
	fftwf_plan forward_ps;
	float* input_fft_buffer_g;
	float* output_fft_buffer_g;
	fftwf_plan forward_g;
	fftwf_plan backward_g;

	//Arrays and variables for getting bins info
	float* fft_p2; //power spectrum
	float* fft_magnitude; //magnitude spectrum
	float* fft_phase; //phase spectrum

	//noise related
	float* noise_thresholds_p2; //captured noise profile power spectrum
	float* noise_thresholds_scaled; //captured noise profile power spectrum scaled by oversustraction
	float* spectral_envelope_values; //spectral envelope
	FFTPeak* noise_spectral_peaks; //Tonal noises of the noise profile
	int* peak_pos; //Whether theres a peak or not in current postiion of the spectrum
	int peak_count;	//counts the amount of peaks found
	bool peaks_detected; //flag to execute peak detection

	//smoothing related
	float* smoothed_spectrum; //power spectrum to be smoothed
	float* smoothed_spectrum_prev; //previous frame smoothed power spectrum for envelopes

	//Transient preservation related
	float* transient_preserv_prev; //previous frame smoothed power spectrum for envelopes

	//Reduction gains
	float* Gk_spectral_subtaction; //gain to be applied
	float* Gk_spectral_gates; //gain to be applied
	float* Gk; //definitive gain

	//Ensemble related
	float* residual_spectrum;
	float* denoised_spectrum;
	float* final_spectrum;

	//Loizou algorithm
	float* auto_thresholds; //Reference threshold for louizou algorithm
	float* prev_noise_thresholds;
	float* s_pow_spec;
	float* prev_s_pow_spec;
	float* p_min;
	float* prev_p_min;
	float* speech_p_p;
	float* prev_speech_p_p;

	//masking
	float* bark_z;
	float* absolute_thresholds; //absolute threshold of hearing
	float* SSF;
	float* spl_reference_values;
	float* unity_gain_bark_spectrum;
	float* spreaded_unity_gain_bark_spectrum;
	float* alpha_masking;
	float* beta_masking;
	float* input_fft_buffer_at;
	float* output_fft_buffer_at;
	fftwf_plan forward_at;

	// clock_t start, end;
	// double cpu_time_used;

	//LV2 state URID (Save and restore noise profile)
	LV2_URID_Map* map;
	LV2_URID atom_Vector;
	LV2_URID atom_Int;
	LV2_URID atom_Float;
	LV2_URID prop_fftsize;
	LV2_URID prop_nwindow;
	LV2_URID prop_nrepelFFTp2;
} Nrepel;

static LV2_Handle
instantiate(const LV2_Descriptor* descriptor, double rate, const char* bundle_path,
						const LV2_Feature* const* features)
{
	//Actual struct declaration
	Nrepel* nrepel = (Nrepel*)malloc(sizeof(Nrepel));

	//Retrieve the URID map callback, and needed URIDs
	for (int i=0; features[i]; ++i)
	{
		if (!strcmp(features[i]->URI, LV2_URID__map))
		{
			nrepel->map = (LV2_URID_Map*)features[i]->data;
		}
	}
	if (!nrepel->map)
	{
		//bail out: host does not support urid:map
		free(nrepel);
		return NULL;
	}

	//For lv2 state (noise profile saving)
	nrepel->atom_Vector = nrepel->map->map(nrepel->map->handle, LV2_ATOM__Vector);
	nrepel->atom_Int = nrepel->map->map(nrepel->map->handle, LV2_ATOM__Int);
	nrepel->atom_Float = nrepel->map->map(nrepel->map->handle, LV2_ATOM__Float);
	nrepel->prop_fftsize = nrepel->map->map(nrepel->map->handle, NREPEL_URI "#fftsize");
	nrepel->prop_nwindow = nrepel->map->map(nrepel->map->handle, NREPEL_URI "#nwindow");
	nrepel->prop_nrepelFFTp2 = nrepel->map->map(nrepel->map->handle, NREPEL_URI "#FFTp2");

	//Sampling related
	nrepel->samp_rate = (float)rate;

	//FFT related
	//nrepel->fft_size = next_pow_two(FRAME_SIZE/1000.f * nrepel->samp_rate);
	nrepel->fft_size = FFT_SIZE;
	nrepel->fft_size_2 = nrepel->fft_size/2;
	nrepel->input_fft_buffer = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->output_fft_buffer = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->forward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer, nrepel->output_fft_buffer, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->backward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->output_fft_buffer, nrepel->input_fft_buffer, FFTW_HC2R, FFTW_ESTIMATE);

	//STFT window related
	nrepel->window_option_input = INPUT_WINDOW;
	nrepel->window_option_output = OUTPUT_WINDOW;
	//nrepel->frame_size = nearest_odd(FRAME_SIZE/1000.f * nrepel->samp_rate);//correct phase response
	nrepel->frame_size = FRAME_SIZE;
	nrepel->frame_fft_diff = nrepel->fft_size - nrepel->frame_size;
	nrepel->input_window = (float*)calloc(nrepel->frame_size, sizeof(float));
	nrepel->output_window = (float*)calloc(nrepel->frame_size, sizeof(float));

	//buffers for OLA
	nrepel->in_fifo = (float*)calloc(nrepel->frame_size, sizeof(float));
	nrepel->out_fifo = (float*)calloc(nrepel->frame_size, sizeof(float));
	nrepel->output_accum = (float*)calloc(nrepel->frame_size*2, sizeof(float));
	nrepel->overlap_factor = OVERLAP_FACTOR;
	nrepel->hop = nrepel->frame_size/nrepel->overlap_factor;
	nrepel->input_latency = nrepel->frame_size - nrepel->hop;
	nrepel->read_ptr = nrepel->input_latency; //the initial position because we are that many samples ahead

	//soft bypass
	nrepel->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f  / nrepel->samp_rate));
	nrepel->wet_dry = 0.f;

	//Arrays for getting bins info
	nrepel->fft_p2 = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->fft_magnitude = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->fft_phase = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));

	//noise threshold related
	nrepel->noise_thresholds_p2 = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->noise_thresholds_scaled = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->spectral_envelope_values = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->noise_spectral_peaks = (FFTPeak*)calloc((SP_MAX_NUM), sizeof(FFTPeak));
	nrepel->peak_pos = (int*)calloc((nrepel->fft_size_2+1), sizeof(int));
	nrepel->window_count = 0.f;
	nrepel->peak_count = 0.f;
	nrepel->peak_detection = false;
	nrepel->noise_thresholds_availables = false;

	//noise adaptive estimation related
	nrepel->auto_thresholds = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->prev_noise_thresholds = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->s_pow_spec = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->prev_s_pow_spec = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->p_min = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->prev_p_min = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->speech_p_p = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->prev_speech_p_p = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));

	//smoothing related
	nrepel->smoothed_spectrum = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->smoothed_spectrum_prev = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->prev_beta = 0.f;

	//transient preservation
	nrepel->transient_preserv_prev = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));


	//masking related
	nrepel->bark_z = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->absolute_thresholds = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->unity_gain_bark_spectrum = (float*)calloc(N_BARK_BANDS, sizeof(float));
	nrepel->spreaded_unity_gain_bark_spectrum = (float*)calloc(N_BARK_BANDS, sizeof(float));
	nrepel->spl_reference_values = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->alpha_masking = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->beta_masking = (float*)calloc((nrepel->fft_size_2+1), sizeof(float));
	nrepel->SSF = (float*)calloc((N_BARK_BANDS*N_BARK_BANDS), sizeof(float));
	nrepel->input_fft_buffer_at = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->output_fft_buffer_at = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->forward_at = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer_at, nrepel->output_fft_buffer_at, FFTW_R2HC, FFTW_ESTIMATE);

	//reduction gains related
	nrepel->Gk_spectral_subtaction = (float*)calloc((nrepel->fft_size), sizeof(float));
	nrepel->Gk_spectral_gates = (float*)calloc((nrepel->fft_size), sizeof(float));
	nrepel->Gk = (float*)calloc((nrepel->fft_size), sizeof(float));

	//Postfilter related
	nrepel->input_fft_buffer_ps = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->output_fft_buffer_ps = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->forward_ps = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer_ps, nrepel->output_fft_buffer_ps, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->input_fft_buffer_g = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->output_fft_buffer_g = (float*)calloc(nrepel->fft_size, sizeof(float));
	nrepel->forward_g = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer_g, nrepel->output_fft_buffer_g, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->backward_g = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->output_fft_buffer_g, nrepel->input_fft_buffer_g, FFTW_HC2R, FFTW_ESTIMATE);

	//final ensemble related
	nrepel->residual_spectrum = (float*)calloc((nrepel->fft_size), sizeof(float));
	nrepel->denoised_spectrum = (float*)calloc((nrepel->fft_size), sizeof(float));
	nrepel->final_spectrum = (float*)calloc((nrepel->fft_size), sizeof(float));

	//Window combination initialization (pre processing window post processing window)
	fft_pre_and_post_window(nrepel->input_window, nrepel->output_window,
													nrepel->frame_size, nrepel->window_option_input,
													nrepel->window_option_output,	&nrepel->overlap_scale_factor);

	//Set initial gain as unity for the positive part
	initialize_array(nrepel->Gk_spectral_subtaction,1.f,nrepel->fft_size);
	initialize_array(nrepel->Gk_spectral_gates,1.f,nrepel->fft_size);
	initialize_array(nrepel->Gk,1.f,nrepel->fft_size);

	//Compute adaptive initial thresholds
	compute_auto_thresholds(nrepel->auto_thresholds, nrepel->fft_size, nrepel->fft_size_2,
													nrepel->samp_rate);

	//MASKING initializations
	compute_bark_mapping(nrepel->bark_z, nrepel->fft_size_2, nrepel->samp_rate);
	compute_absolute_thresholds(nrepel->absolute_thresholds, nrepel->fft_size_2,
															nrepel->samp_rate);
	spl_reference(nrepel->spl_reference_values, nrepel->fft_size_2, nrepel->samp_rate,
								nrepel->input_fft_buffer_at, nrepel->output_fft_buffer_at,
								&nrepel->forward_at);
	compute_SSF(nrepel->SSF);

	//Initializing unity gain values for offset normalization
	initialize_array(nrepel->unity_gain_bark_spectrum,1.f,N_BARK_BANDS);
	//Convolve unitary energy bark spectrum with SSF
  convolve_with_SSF(nrepel->SSF, nrepel->unity_gain_bark_spectrum,
										nrepel->spreaded_unity_gain_bark_spectrum);

	initialize_array(nrepel->alpha_masking,1.f,nrepel->fft_size_2+1);
	initialize_array(nrepel->beta_masking,0.f,nrepel->fft_size_2+1);

	return (LV2_Handle)nrepel;
}

static void
connect_port(LV2_Handle instance, uint32_t port, void* data)
{
	Nrepel* nrepel = (Nrepel*)instance;

	switch ((PortIndex)port)
	{
		case NREPEL_AMOUNT:
		nrepel->amount_of_reduction = (float*)data;
		break;
		case NREPEL_NOFFSET:
		nrepel->noise_thresholds_offset = (float*)data;
		break;
		case NREPEL_RELEASE:
		nrepel->release = (float*)data;
		break;
		case NREPEL_MASKING:
		nrepel->masking = (float*)data;
		break;
		case NREPEL_A_CONTROL:
		nrepel->artifact_control = (float*)data;
		break;
		case NREPEL_WHITENING:
		nrepel->whitening_factor_pc = (float*)data;
		break;
		case NREPEL_N_LEARN:
		nrepel->noise_learn_state = (float*)data;
		break;
		case NREPEL_N_ADAPTIVE:
		nrepel->adaptive_state = (float*)data;
		break;
		case NREPEL_RESIDUAL_LISTEN:
		nrepel->residual_listen = (float*)data;
		break;
		case NREPEL_RESET:
		nrepel->reset_profile = (float*)data;
		break;
		case NREPEL_ENABLE:
		nrepel->enable = (float*)data;
		break;
		case NREPEL_LATENCY:
		nrepel->report_latency = (float*)data;
		break;
		case NREPEL_INPUT:
		nrepel->input = (const float*)data;
		break;
		case NREPEL_OUTPUT:
		nrepel->output = (float*)data;
		break;
	}
}

static void
activate(LV2_Handle instance)
{
}

static void
run(LV2_Handle instance, uint32_t n_samples)
{
	Nrepel* nrepel = (Nrepel*)instance;

	// //Time execution measurement
	// nrepel->start = clock();
	// //--------------

	//handy variables
	int k;
	unsigned int pos;

	//Inform latency at run call
	*(nrepel->report_latency) = (float) nrepel->input_latency;

	//Reset button state (if on)
	if (*(nrepel->reset_profile) == 1.f)
	{
		initialize_array(nrepel->noise_thresholds_p2,0.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->noise_thresholds_scaled,0.f,nrepel->fft_size_2+1);
		nrepel->window_count = 0.f;
		nrepel->peak_count = 0.f;
		nrepel->noise_thresholds_availables = false;
		nrepel->peak_detection = false;

		initialize_array(nrepel->Gk_spectral_subtaction,1.f,nrepel->fft_size);
		initialize_array(nrepel->Gk_spectral_gates,1.f,nrepel->fft_size);
		initialize_array(nrepel->Gk,1.f,nrepel->fft_size);

		initialize_array(nrepel->prev_noise_thresholds,0.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->s_pow_spec,0.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->prev_s_pow_spec,0.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->p_min,0.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->prev_p_min,0.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->speech_p_p,0.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->prev_speech_p_p,0.f,nrepel->fft_size_2+1);

		initialize_array(nrepel->alpha_masking,1.f,nrepel->fft_size_2+1);
		initialize_array(nrepel->beta_masking,0.f,nrepel->fft_size_2+1);
	}

	//Softbypass targets in case of disabled or enabled
	if(*(nrepel->enable) == 0.f)
	{ //if disabled
		nrepel->wet_dry_target = 0.f;
	}
	else
	{ //if enabled
		nrepel->wet_dry_target = 1.f;
	}
	//Interpolate parameters over time softly to bypass without clicks or pops
	nrepel->wet_dry += nrepel->tau * (nrepel->wet_dry_target - nrepel->wet_dry) + FLT_MIN;

	//Parameters values
	/*exponential decay coefficients for envelopes and adaptive noise profiling
		These must take into account the hop size as explained in the following paper
		FFT-BASED DYNAMIC RANGE COMPRESSION*/
	nrepel->release_coeff = expf(-1000.f/(((*(nrepel->release)) * nrepel->samp_rate)/ nrepel->hop) );
	nrepel->amount_of_reduction_linear = from_dB(-1.f * *(nrepel->amount_of_reduction));
	nrepel->thresholds_offset_linear = from_dB(*(nrepel->noise_thresholds_offset));
	nrepel->whitening_factor = *(nrepel->whitening_factor_pc)/100.f;

	//main loop for processing
	for (pos = 0; pos < n_samples; pos++)
	{
		//Store samples int the input buffer
		nrepel->in_fifo[nrepel->read_ptr] = nrepel->input[pos];
		//Output samples in the output buffer (even zeros introduced by latency)
		nrepel->output[pos] = nrepel->out_fifo[nrepel->read_ptr - nrepel->input_latency];
		//Now move the read pointer
		nrepel->read_ptr++;

		//Once the buffer is full we can do stuff
		if (nrepel->read_ptr >= nrepel->frame_size)
		{
			//Reset the input buffer position
			nrepel->read_ptr = nrepel->input_latency;

			//----------STFT Analysis------------

			int start_pos = nrepel->frame_fft_diff/2;

			//Adding and windowing the frame input values in the center (zero-phasing)
			for (k = 0; k < nrepel->frame_size; k++)
			{
				nrepel->input_fft_buffer[k] = nrepel->in_fifo[start_pos + k] * nrepel->input_window[k];
			}

			//----------FFT Analysis------------

			//Do transform
			fftwf_execute(nrepel->forward);

			//-----------GET INFO FROM BINS--------------

			get_info_from_bins(nrepel->fft_p2, nrepel->fft_magnitude, nrepel->fft_phase,
												 nrepel->fft_size_2, nrepel->fft_size,
												 nrepel->output_fft_buffer);

			/////////////////////SPECTRAL PROCESSING//////////////////////////

			/*This section countains the specific noise reduction processing blocks
				but it could be replaced with any spectral processing (I'm looking at you future tinkerer)
				Parameters for the STFT transform can be changed at the top of this file
			*/

			//If the spectrum is not silence
			if(!is_empty(nrepel->fft_p2, nrepel->fft_size_2))
			{
				//If adaptive noise is selected the noise is adapted in time
				if(*(nrepel->adaptive_state) == 1.f)
				{
					//This has to be revised(issue 8 on github)
					adapt_noise(nrepel->fft_p2, nrepel->fft_size_2, nrepel->noise_thresholds_p2,
											nrepel->auto_thresholds, nrepel->prev_noise_thresholds,
											nrepel->s_pow_spec, nrepel->prev_s_pow_spec, nrepel->p_min,
											nrepel->prev_p_min, nrepel->speech_p_p, nrepel->prev_speech_p_p);

					nrepel->noise_thresholds_availables = true;
				}

				/*If selected estimate noise spectrum is based on selected portion of signal
				 *do not process the signal
				 */
				if(*(nrepel->noise_learn_state) == 1.f)
				{ //MANUAL
					get_noise_statistics(nrepel->fft_p2, nrepel->fft_size_2,
															 nrepel->noise_thresholds_p2, &nrepel->window_count);

					nrepel->noise_thresholds_availables = true;
				}
				else
				{
					//If there is a noise profile reduce noise
					if (nrepel->noise_thresholds_availables == true)
					{
						//Detector smoothing and oversustraction
						preprocessing(nrepel->thresholds_offset_linear, nrepel->fft_p2,
													nrepel->noise_thresholds_p2, nrepel->noise_thresholds_scaled,
													nrepel->smoothed_spectrum, nrepel->smoothed_spectrum_prev,
													nrepel->fft_size_2, &nrepel->prev_beta, nrepel->bark_z,
													nrepel->absolute_thresholds, nrepel->SSF,
													nrepel->release_coeff,
													nrepel->spreaded_unity_gain_bark_spectrum,
													nrepel->spl_reference_values, nrepel->alpha_masking,
													nrepel->beta_masking, *(nrepel->masking));

						//Supression rule
						spectral_gain(nrepel->fft_p2, nrepel->noise_thresholds_p2, nrepel->noise_thresholds_scaled, nrepel->smoothed_spectrum,
													nrepel->fft_size_2, *(nrepel->adaptive_state), nrepel->Gk, nrepel->Gk_spectral_subtaction, nrepel->Gk_spectral_gates,
													*(nrepel->artifact_control), nrepel->transient_preserv_prev);

						//apply gains
						denoised_calulation(nrepel->fft_size, nrepel->output_fft_buffer,
																nrepel->denoised_spectrum, nrepel->Gk);

						//residual signal
						residual_calulation(nrepel->fft_size, nrepel->output_fft_buffer,
																nrepel->residual_spectrum, nrepel->denoised_spectrum,
																nrepel->whitening_factor);

						//Ensemble the final spectrum using residual and denoised
						final_spectrum_ensemble(nrepel->fft_size,nrepel->final_spectrum,
																		nrepel->residual_spectrum,
																		nrepel->denoised_spectrum,
																		nrepel->amount_of_reduction_linear,
																		*(nrepel->residual_listen));

						soft_bypass(nrepel->final_spectrum, nrepel->output_fft_buffer,
												nrepel->wet_dry, nrepel->fft_size);
					}
				}
			}

			///////////////////////////////////////////////////////////

			//----------STFT Synthesis------------

			//------------FFT Synthesis-------------

			//Do inverse transform
			fftwf_execute(nrepel->backward);

			//Normalizing value
			for (k = 0; k < nrepel->fft_size; k++)
			{
				nrepel->input_fft_buffer[k] = nrepel->input_fft_buffer[k]/nrepel->fft_size;
			}

			//------------OVERLAPADD-------------

			//Windowing scaling and accumulation
			for(k = 0; k < nrepel->frame_size; k++)
			{
				nrepel->output_accum[k] += (nrepel->output_window[k]*nrepel->input_fft_buffer[start_pos + k])/(nrepel->overlap_scale_factor * nrepel->overlap_factor);
			}

			//Output samples up to the hop size
			for (k = 0; k < nrepel->hop; k++)
			{
				nrepel->out_fifo[k] = nrepel->output_accum[k];
			}

			//shift FFT accumulator the hop size
			memmove(nrepel->output_accum, nrepel->output_accum + nrepel->hop,
							nrepel->frame_size*sizeof(float));

			//move input FIFO
			for (k = 0; k < nrepel->input_latency; k++)
			{
				nrepel->in_fifo[k] = nrepel->in_fifo[k+nrepel->hop];
			}
			//-------------------------------
		}//if
	}//main loop

	// //Time measurement
	// nrepel->end = clock();
	// nrepel->cpu_time_used = ((double) (nrepel->end - nrepel->start)) / CLOCKS_PER_SEC;
	//
	// //To string
	// char buffer[50];
	// sprintf(buffer,"%lf",nrepel->cpu_time_used);
	// strcat(buffer,"\n");
	//
	// //Saving results to a file
	// FILE *fp;
	//
	// fp = fopen("resuts.txt", "a");
	// fputs(buffer, fp);
	// fclose(fp);
}

static void
deactivate(LV2_Handle instance)
{
}

static void
cleanup(LV2_Handle instance)
{
	free(instance);
}

static LV2_State_Status
savestate(LV2_Handle instance, LV2_State_Store_Function store, LV2_State_Handle handle,
					uint32_t flags, const LV2_Feature* const* features)
{
	Nrepel* nrepel = (Nrepel*)instance;

	FFTVector* vector = (FFTVector*)malloc(sizeof(FFTVector));

	vector->child_type = nrepel->atom_Float;
	vector->child_size = sizeof(float);

	store(handle, nrepel->prop_fftsize, &nrepel->fft_size_2, sizeof(int), nrepel->atom_Int,
				LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	store(handle, nrepel->prop_nwindow, &nrepel->window_count, sizeof(float),
				nrepel->atom_Float, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	memcpy(vector->array, nrepel->noise_thresholds_p2, sizeof(vector->array));

	store(handle, nrepel->prop_nrepelFFTp2, (void*) vector, sizeof(FFTVector),
				nrepel->atom_Vector, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

  return LV2_STATE_SUCCESS;
}

static LV2_State_Status
restorestate(LV2_Handle instance, LV2_State_Retrieve_Function retrieve,
						 LV2_State_Handle handle, uint32_t flags,
						 const LV2_Feature* const* features)
{
	Nrepel* nrepel = (Nrepel*)instance;
	size_t   size;
	uint32_t type;
	uint32_t valflags;

	const int32_t* fftsize = retrieve(handle, nrepel->prop_fftsize, &size, &type, &valflags);
	if (!fftsize || type != nrepel->atom_Int || *fftsize != nrepel->fft_size_2){
		return LV2_STATE_ERR_NO_PROPERTY;
	}

	const void* vecFFTp2 = retrieve(handle, nrepel->prop_nrepelFFTp2, &size, &type, &valflags);
	if ( !vecFFTp2 || size != sizeof(FFTVector) || type != nrepel->atom_Vector){
		return LV2_STATE_ERR_NO_PROPERTY;
	}

	//Deactivate any denoising before loading any noise profile
	nrepel->noise_thresholds_availables = false;

	//Copy to local variables
	memcpy(nrepel->noise_thresholds_p2, (float*) LV2_ATOM_BODY(vecFFTp2), (nrepel->fft_size_2+1)*sizeof(float));

	const float* wincount = retrieve(handle, nrepel->prop_nwindow, &size, &type, &valflags);
	if (fftsize && type == nrepel->atom_Float) {
		nrepel->window_count = *wincount;
	}

	//Reactivate denoising with restored profile
	nrepel->noise_thresholds_availables = true;

	return LV2_STATE_SUCCESS;
}

static const void*
extension_data(const char* uri)
{
	static const LV2_State_Interface  state  = { savestate, restorestate };
	if (!strcmp(uri, LV2_STATE__interface))
	{
		return &state;
	}
	return NULL;
}

static const
LV2_Descriptor descriptor =
{
	NREPEL_URI,
	instantiate,
	connect_port,
	activate,
	run,
	deactivate,
	cleanup,
	extension_data
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor*
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
