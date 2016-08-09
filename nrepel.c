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
#include <stdio.h>
#include <string.h>
#include <fftw3.h>


#include "denoise_gain.c"

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"

//Noise capture states
#define OFF_STATE 0
#define MANUAL_CAPTURE_ON_STATE 1
#define ADAPTIVE_CAPTURE_STATE 2

//STFT default values
#define DEFAULT_FFT_SIZE 2048 //max should be 8192 otherwise is too expensive
#define DEFAULT_WINDOW_TYPE 2 //0 Hann 1 Hamming 2 Blackman
#define DEFAULT_OVERLAP_FACTOR 4 //2 is 50% and 4 is 75% overlap

//Denoise related options
// #define NOISE_STAT_CHOISE 1 //0 max 1 geometric mean 2 average
// #define DENOISE_METHOD 2 //0 Wiener 1 Power Substraction 2 EM with CM
#define ALPHA 0.98 //For EM denoising
#define WA 0.05 //For spectral whitening strenght 0-1

///---------------------------------------------------------------------

//LV2 CODE

typedef enum {
	NREPEL_CAPTURE = 0,
	NREPEL_N_THRESH = 1,
	NREPEL_SMOOTHING = 2,
	NREPEL_METHOD = 3,
	NREPEL_AMOUNT = 4,
	NREPEL_OVERRED = 5,
	NREPEL_LATENCY = 6,
	NREPEL_WHITENING = 7,
	NREPEL_RESET = 8,
	NREPEL_NOISE_LISTEN = 9,
	NREPEL_INPUT = 10,
	NREPEL_OUTPUT = 11,
} PortIndex;

typedef struct {
	const float* input; //input of samples from host (changing size)
	float* output; //output of samples to host (changing size)
	float samp_rate; // Sample rate received from the host

	//Parameters for the algorithm (user input)
	float* capt_state; // Capture Noise state (Manual-Off-Auto)
	float* amount_reduc; // Amount of noise to reduce in dB
	float* over_reduc; // Amount of noise to reduce in dB
	float* report_latency; // Latency necessary
	float* reset_print; // Latency necessary
	float* noise_listen; //For noise only listening
	float* noise_stat_choise; //Choise of statistic to use for noise spectrum
	float* denoise_method; //Choise of denoise method
	float* residue_whitening; //Whitening of the residual spectrum
	float* noise_thresh; //detection threshold for louizou estimation

	//Parameters for the STFT
	int fft_size; //FFTW input size
	int fft_size_2; //FFTW half input size
	int window_type; // Window type for the STFT
	int overlap_factor; //oversampling factor for overlap calculations
	int hop; // Hop size for the STFT
	float* window; // Window values

	//Temporary buffers for processing and outputting
	int input_latency;
	float* in_fifo; //internal input buffer
	float* out_fifo; //internal output buffer
	float* output_accum; //FFT output accumulator
	int read_ptr; //buffers read pointer

	//FFTW related variables
	float* input_fft_buffer;
	float* output_fft_buffer;
	fftwf_plan forward;
	fftwf_plan backward;

	float real_p,real_n,mag,p2;
	float* fft_magnitude;//magnitude
	float* fft_p2;//power
	float* fft_p2_prev;//power previous frame

	//Store variables
	float max_float;
	float* noise_print_min; // The min noise spectrum computed by the captured signal
	float* noise_print_max; // The max noise spectrum computed by the captured signal
	float* noise_spectrum;
	float* residual_spectrum;
	float wa;
	float* tappering_filter;
	float* whitening_spectrum;


	float* Gk; //gain to be applied
	float* Gk_prev; //last gain applied
	float* gain_prev; //previously gain computed
	float alpha;//for EM
	int prev_frame;//for EM


	float* prev_noise;
	float* s_pow_spec;
	float* prev_s_pow_spec;
	float* p_min;
	float* prev_p_min;
	float* speech_p_p;
	float* prev_speech_p_p;

	float sum_log_p;
	float sum_p;
	float kinv;
	float SFM;
	float tonality_factor;
	float* masked;
	float** jg_upper;
	float** jg_lower;
	float* bark_z;

} Nrepel;

static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
						double                    rate,
						const char*               bundle_path,
						const LV2_Feature* const* features) {
	//Actual struct declaration
	Nrepel* nrepel = (Nrepel*)malloc(sizeof(Nrepel));

	//Initialize variables
	nrepel->samp_rate = rate;
	nrepel->fft_size = DEFAULT_FFT_SIZE;
	nrepel->window_type = DEFAULT_WINDOW_TYPE;
	nrepel->overlap_factor = DEFAULT_OVERLAP_FACTOR;
	nrepel->max_float = FLT_MAX;
	nrepel->alpha = ALPHA;
	nrepel->wa = WA;
	nrepel->prev_frame = 0;

	nrepel->fft_size_2 = nrepel->fft_size/2;
	nrepel->hop = nrepel->fft_size/nrepel->overlap_factor;
	nrepel->input_latency = nrepel->fft_size - nrepel->hop;
	nrepel->read_ptr = nrepel->input_latency; //the initial position because we are that many samples ahead
	nrepel->kinv	= 1.f/(float)nrepel->fft_size_2;

	nrepel->in_fifo = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->out_fifo = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_accum = (float*)calloc(nrepel->fft_size,sizeof(float));

	nrepel->window = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->input_fft_buffer = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_fft_buffer = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->forward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer, nrepel->output_fft_buffer, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->backward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->output_fft_buffer, nrepel->input_fft_buffer, FFTW_HC2R, FFTW_ESTIMATE);

	nrepel->fft_magnitude = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_p2 = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_p2_prev = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	nrepel->noise_print_min = (float*)malloc((nrepel->fft_size_2+1)*sizeof(float));
	nrepel->noise_print_max = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->noise_spectrum = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->residual_spectrum = (float*)calloc((nrepel->fft_size),sizeof(float));

	nrepel->Gk = (float*)malloc((nrepel->fft_size_2+1)*sizeof(float));
	nrepel->Gk_prev = (float*)malloc((nrepel->fft_size_2+1)*sizeof(float));
	nrepel->gain_prev = (float*)malloc((nrepel->fft_size_2+1)*sizeof(float));
	nrepel->whitening_spectrum = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->tappering_filter = (float*)calloc(nrepel->fft_size_2+1,sizeof(float));

	nrepel->prev_noise = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->s_pow_spec = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_s_pow_spec = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->p_min = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_p_min = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->speech_p_p = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_speech_p_p = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	nrepel->masked = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->jg_upper = (float**)malloc(sizeof(float)*(nrepel->fft_size));
	nrepel->jg_lower = (float**)malloc(sizeof(float)*(nrepel->fft_size));
	nrepel->bark_z = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	for(int k = 0; k < nrepel->fft_size; k++){
		nrepel->jg_upper[k] = (float*)calloc(11,sizeof(float));
		nrepel->jg_lower[k] = (float*)calloc(11,sizeof(float));
	}

	//Here we initialize arrays with intended default values
	memset(nrepel->noise_print_min, nrepel->max_float, (nrepel->fft_size_2+1)*sizeof(float));
	memset(nrepel->Gk, 1, (nrepel->fft_size_2+1)*sizeof(float));
	memset(nrepel->Gk_prev, 1, (nrepel->fft_size_2+1)*sizeof(float));
	memset(nrepel->gain_prev, 1, (nrepel->fft_size_2+1)*sizeof(float));

	fft_window(nrepel->window,nrepel->fft_size,nrepel->window_type); //STFT window
	tappering_filter_calc(nrepel->tappering_filter,(nrepel->fft_size_2+1),WA); //Tappering window

	//MASKING THRESHOLDS (this should be precomputed and in a file to be more efficient)
	compute_bark_z(nrepel->bark_z,nrepel->fft_size_2,nrepel->samp_rate);

	return (LV2_Handle)nrepel;
}



static void
connect_port(LV2_Handle instance,
						 uint32_t   port,
						 void*      data) {
	Nrepel* nrepel = (Nrepel*)instance;

	switch ((PortIndex)port) {
		case NREPEL_CAPTURE:
		nrepel->capt_state = (float*)data;
		break;
		case NREPEL_N_THRESH:
		nrepel->noise_thresh = (float*)data;
		break;
		case NREPEL_SMOOTHING:
		nrepel->noise_stat_choise = (float*)data;
		break;
		case NREPEL_METHOD:
		nrepel->denoise_method = (float*)data;
		break;
		case NREPEL_AMOUNT:
		nrepel->amount_reduc = (float*)data;
		break;
		case NREPEL_OVERRED:
		nrepel->over_reduc = (float*)data;
		break;
		case NREPEL_LATENCY:
		nrepel->report_latency = (float*)data;
		break;
		case NREPEL_WHITENING:
		nrepel->residue_whitening = (float*)data;
		break;
		case NREPEL_RESET:
		nrepel->reset_print = (float*)data;
		break;
		case NREPEL_NOISE_LISTEN:
		nrepel->noise_listen = (float*)data;
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
run(LV2_Handle instance, uint32_t n_samples) {
	Nrepel* nrepel = (Nrepel*)instance;

	//handy variables
	int k;
	unsigned int pos;

	//Inform latency at run call
	*(nrepel->report_latency) = (float) nrepel->input_latency;

	//Reset reset button state
	if (*(nrepel->reset_print) == 1.f) {
		memset(nrepel->noise_print_min, nrepel->max_float, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->noise_print_max, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->noise_spectrum, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->residual_spectrum, 0, (nrepel->fft_size)*sizeof(float));
		memset(nrepel->Gk, 1, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->Gk_prev, 1, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->gain_prev, 1, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_noise, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->s_pow_spec, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_s_pow_spec, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->p_min, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_p_min, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->speech_p_p, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_speech_p_p, 0, (nrepel->fft_size_2+1)*sizeof(float));
		nrepel->prev_frame = 0;
	}

	//main loop for processing
	for (pos = 0; pos < n_samples; pos++){
		//Store samples int the input buffer
		nrepel->in_fifo[nrepel->read_ptr] = nrepel->input[pos];
		//Output samples in the output buffer (even zeros introduced by latency)
		nrepel->output[pos] = nrepel->out_fifo[nrepel->read_ptr - nrepel->input_latency];
		//Now move the read pointer
		nrepel->read_ptr++;

		//Once the buffer is full we can do stuff
		if (nrepel->read_ptr >= nrepel->fft_size){
			//Reset the input buffer position
			nrepel->read_ptr = nrepel->input_latency;

			//Apply windowing (Could be any window type)
			for (k = 0; k < nrepel->fft_size; k++){
				nrepel->input_fft_buffer[k] = nrepel->in_fifo[k] * nrepel->window[k];
			}

			//----------FFT Analysis------------
			//Do transform
			fftwf_execute(nrepel->forward);

			//-----------GET INFO FROM BINS--------------

			//Get the positive spectrum and compute the power spectrum or magnitude
			for (k = 0; k <= nrepel->fft_size_2; k++){
				nrepel->real_p = nrepel->output_fft_buffer[k];
				nrepel->real_n = nrepel->output_fft_buffer[nrepel->fft_size-k];

				//Get mag and/or power
				if(k < nrepel->fft_size){
					nrepel->p2 = nrepel->real_p*nrepel->real_p + nrepel->real_n*nrepel->real_n;
					//nrepel->mag = sqrtf(nrepel->p2);
				} else {
					//Nyquist
					nrepel->p2 = nrepel->real_p*nrepel->real_p;
					//nrepel->mag = nrepel->real_p;
				}

				//Store values in magnitude and power arrays
				nrepel->fft_p2[k] = sanitize_denormal(nrepel->p2);
				//nrepel->fft_magnitude[k] = sanitize_denormal(nrepel->mag);

				//FOR MASKING THRESHOLDS
				if (*(nrepel->denoise_method) == 4.f) {
					nrepel->sum_log_p += log10f(nrepel->fft_p2[k]);
			    nrepel->sum_p += nrepel->fft_p2[k];
				}
			}

			//FOR MASKING THRESHOLDS
			if (*(nrepel->denoise_method) == 4.f) {
				nrepel->SFM = 10.f*( nrepel->kinv*nrepel->sum_log_p - log10f(nrepel->sum_p*nrepel->kinv));
				nrepel->tonality_factor = MIN(nrepel->SFM/-60.f, 1.f);
				compute_johnston_gain(nrepel->bark_z,nrepel->jg_upper,nrepel->jg_lower,nrepel->fft_size_2,nrepel->tonality_factor);
			}

			//More info could be useful like onsets for example to modify the amount of reduction

			//------------Processing---------------

			switch ((int) *(nrepel->capt_state) ) {
				case MANUAL_CAPTURE_ON_STATE:
					//If selected estimate noise spectrum based on selected portion of signal
					estimate_noise_spectrum_manual(nrepel->fft_p2,
																	nrepel->fft_size_2,
																	nrepel->noise_spectrum,
																	nrepel->noise_print_min,
					 											  nrepel->noise_print_max,
																	*(nrepel->noise_stat_choise),
																	nrepel->wa,
		                              nrepel->whitening_spectrum);
					break;
				case ADAPTIVE_CAPTURE_STATE:
					//if slected auto estimate noise spectrum and apply denoising
					estimate_noise_spectrum_adaptive(nrepel->fft_p2,
																					 nrepel->fft_size_2,
																					 nrepel->noise_spectrum,
																					 from_dB(*(nrepel->noise_thresh)),
																					 nrepel->prev_noise,
																					 nrepel->s_pow_spec,
																					 nrepel->prev_s_pow_spec,
																					 nrepel->p_min,
																					 nrepel->prev_p_min,
																					 nrepel->speech_p_p,
																					 nrepel->prev_speech_p_p,
																					 nrepel->wa,
																					 nrepel->whitening_spectrum);
					//Compute denoising gain based on previously computed spectrum
					denoise_gain(*(nrepel->denoise_method),
											 *(nrepel->over_reduc),
											 nrepel->fft_p2,
											 nrepel->fft_p2_prev,
											 nrepel->fft_size_2,
											 nrepel->Gk,
											 nrepel->Gk_prev,
											 nrepel->gain_prev,
											 nrepel->noise_spectrum,
											 nrepel->alpha,
											 &nrepel->prev_frame,
										 	 nrepel->masked,
											 nrepel->jg_upper,
		                   nrepel->jg_lower);

					 //APPLY REDUCTION

 					//Residual signal
 					for (k = 0; k <= nrepel->fft_size_2; k++) {
 					 nrepel->residual_spectrum[k] = nrepel->output_fft_buffer[k] - (nrepel->output_fft_buffer[k] * nrepel->Gk[k]);
 					 if(k < nrepel->fft_size_2)
 					 	nrepel->residual_spectrum[nrepel->fft_size-k] = nrepel->output_fft_buffer[nrepel->fft_size-k] - (nrepel->output_fft_buffer[nrepel->fft_size-k] * nrepel->Gk[k]);
 					}

 					//Residue Whitening and tappering
 					if(*(nrepel->residue_whitening) == 1.f) {
 						for (k = 0; k <= nrepel->fft_size_2; k++) {
 							if(nrepel->residual_spectrum[k] > FLT_MIN) {
 								nrepel->residual_spectrum[k] /= nrepel->whitening_spectrum[k];
 								nrepel->residual_spectrum[k] *= nrepel->tappering_filter[k];//Half hann window tappering in favor of high frequencies
 								if(k < nrepel->fft_size_2) {
 									nrepel->residual_spectrum[nrepel->fft_size-k] /= nrepel->whitening_spectrum[k];
 									nrepel->residual_spectrum[nrepel->fft_size-k] *= nrepel->tappering_filter[k];//Half hann window tappering in favor of high frequencies
 								}
 							}
 						}
 					}

 					if (*(nrepel->noise_listen) == 0.f){
 					 //Apply the computed gain to the signal and Mix residual and processed
 					 for (k = 0; k <= nrepel->fft_size_2; k++) {
 						 nrepel->output_fft_buffer[k] *= nrepel->Gk[k];
 						 nrepel->output_fft_buffer[k] += nrepel->residual_spectrum[k]*(1.f/from_dB(*(nrepel->amount_reduc)));
 						 if(k < nrepel->fft_size_2)
 						 	nrepel->output_fft_buffer[nrepel->fft_size-k] *= nrepel->Gk[k];
 						 	nrepel->output_fft_buffer[nrepel->fft_size-k] += nrepel->residual_spectrum[nrepel->fft_size-k]*(1.f/from_dB(*(nrepel->amount_reduc)));
 					 }
 					} else {
 					 //Output noise only
 					 for (k = 0; k <= nrepel->fft_size_2; k++) {
 						 nrepel->output_fft_buffer[k] = nrepel->residual_spectrum[k];
 						 if(k < nrepel->fft_size_2)
 							 nrepel->output_fft_buffer[nrepel->fft_size-k] = nrepel->residual_spectrum[k];
 					 }
 					}
					break;
				case OFF_STATE:
					//Compute denoising gain based on previously computed spectrum (manual or automatic)
					denoise_gain(*(nrepel->denoise_method),
											 *(nrepel->over_reduc),
											 nrepel->fft_p2,
											 nrepel->fft_p2_prev,
											 nrepel->fft_size_2,
											 nrepel->Gk,
											 nrepel->Gk_prev,
											 nrepel->gain_prev,
											 nrepel->noise_spectrum,
											 nrepel->alpha,
											 &nrepel->prev_frame,
											 nrepel->masked,
											 nrepel->jg_upper,
											 nrepel->jg_lower);

					//APPLY REDUCTION

					//Residual signal
					for (k = 0; k <= nrepel->fft_size_2; k++) {
					 nrepel->residual_spectrum[k] = nrepel->output_fft_buffer[k] - (nrepel->output_fft_buffer[k] * nrepel->Gk[k]);
					 if(k < nrepel->fft_size_2)
					 	nrepel->residual_spectrum[nrepel->fft_size-k] = nrepel->output_fft_buffer[nrepel->fft_size-k] - (nrepel->output_fft_buffer[nrepel->fft_size-k] * nrepel->Gk[k]);
					}

					//Residue Whitening and tappering
					if(*(nrepel->residue_whitening) == 1.f) {
						for (k = 0; k <= nrepel->fft_size_2; k++) {
							if(nrepel->residual_spectrum[k] > FLT_MIN) {
								nrepel->residual_spectrum[k] /= nrepel->whitening_spectrum[k];
								nrepel->residual_spectrum[k] *= nrepel->tappering_filter[k];//Half hann window tappering in favor of high frequencies
								if(k < nrepel->fft_size_2) {
									nrepel->residual_spectrum[nrepel->fft_size-k] /= nrepel->whitening_spectrum[k];
									nrepel->residual_spectrum[nrepel->fft_size-k] *= nrepel->tappering_filter[k];//Half hann window tappering in favor of high frequencies
								}
							}
						}
					}

					if (*(nrepel->noise_listen) == 0.f){
					 //Apply the computed gain to the signal and Mix residual and processed
					 for (k = 0; k <= nrepel->fft_size_2; k++) {
						 nrepel->output_fft_buffer[k] *= nrepel->Gk[k];
						 nrepel->output_fft_buffer[k] += nrepel->residual_spectrum[k]*(1.f/from_dB(*(nrepel->amount_reduc)));
						 if(k < nrepel->fft_size_2)
						 	nrepel->output_fft_buffer[nrepel->fft_size-k] *= nrepel->Gk[k];
						 	nrepel->output_fft_buffer[nrepel->fft_size-k] += nrepel->residual_spectrum[nrepel->fft_size-k]*(1.f/from_dB(*(nrepel->amount_reduc)));
					 }
					} else {
					 //Output noise only
					 for (k = 0; k <= nrepel->fft_size_2; k++) {
						 nrepel->output_fft_buffer[k] = nrepel->residual_spectrum[k];
						 if(k < nrepel->fft_size_2)
							 nrepel->output_fft_buffer[nrepel->fft_size-k] = nrepel->residual_spectrum[k];
					 }
					}
					break;
			}

			//------------FFT Synthesis-------------

			//Do inverse transform
			fftwf_execute(nrepel->backward);

			//Scaling FFT (because is not scaled down when backward plan is executed)
			for(k = 0; k < nrepel->fft_size; k++){
				nrepel->input_fft_buffer[k] = nrepel->input_fft_buffer[k]/nrepel->fft_size;
			}

			//Accumulate (Overlapadd)
			for(k = 0; k < nrepel->fft_size; k++){
				nrepel->output_accum[k] += nrepel->input_fft_buffer[k]*nrepel->hop;
			}

			//Output samples up to the hop size
			for (k = 0; k < nrepel->hop; k++){
				nrepel->out_fifo[k] = nrepel->output_accum[k];
			}

			//shift FFT accumulator the hop size
			memmove(nrepel->output_accum, nrepel->output_accum + nrepel->hop, nrepel->fft_size*sizeof(float));

			//Make sure that the non overlaping section is 0
			for (k = (nrepel->fft_size-nrepel->hop); k < nrepel->fft_size; k++){
				nrepel->output_accum[k] = 0.f;
			}

			//move input FIFO
			for (k = 0; k < nrepel->input_latency; k++){
				nrepel->in_fifo[k] = nrepel->in_fifo[k+nrepel->hop];
			}
		}//if
	}//main loop
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

const void*
extension_data(const char* uri)
{
	return NULL;
}

static const
LV2_Descriptor descriptor = {
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
	switch (index) {
		case 0:
		return &descriptor;
		default:
		return NULL;
	}
}
