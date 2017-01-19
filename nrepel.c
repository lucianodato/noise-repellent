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

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"
#define NS_MY "http://example.org/myplugin/schema#"

//STFT default values (This are standart values)
#define FFT_SIZE 2048 //max should be 8192 otherwise is too expensive
#define WINDOW_COMBINATION 0 //0 HANN-HANN 1 HAMMING-HANN 2 BLACKMAN-HANN
#define OVERLAP_FACTOR 4 //4 is 75% overlap
#define HANN_HANN_SCALING 0.375 //This is for overlapadd scaling
#define HAMMING_HANN_SCALING 0.385 // 1/average(window[i]^2)
#define BLACKMAN_HANN_SCALING 0.335

///---------------------------------------------------------------------

//LV2 CODE

typedef enum {
	NREPEL_CAPTURE = 0,
	NREPEL_N_AUTO = 1,
	NREPEL_AMOUNT = 2,
	NREPEL_STRENGTH = 3,
	NREPEL_SMOOTHING = 4,
	// NREPEL_ATTACK = 5,
	// NREPEL_RELEASE = 6,
	NREPEL_FREQUENCY_SMOOTHING = 5,
	NREPEL_MASKING = 6,
	NREPEL_LATENCY = 7,
	NREPEL_WHITENING = 8,
	NREPEL_RESET = 9,
	NREPEL_NOISE_LISTEN = 10,
	NREPEL_ENABLE = 11,
	NREPEL_INPUT = 12,
	NREPEL_OUTPUT = 13,
} PortIndex;

typedef struct {
	const float* input; //input of samples from host (changing size)
	float* output; //output of samples to host (changing size)
	float samp_rate; // Sample rate received from the host

	//Parameters for the algorithm (user input)
	float* capture_state; // Capture Noise state (Manual-Off-Auto)
	float* amount_of_reduction; // Amount of noise to reduce in dB
	float* reduction_strenght; // Amount of noise to reduce in dB
	float* report_latency; // Latency necessary
	float* reset_print; // Latency necessary
	float* noise_listen; //For noise only listening
	float* residual_whitening; //Whitening of the residual spectrum
	float* time_smoothing; //constant that set the time smoothing coefficient
	float* auto_state; //autocapture switch
	float* frequency_smoothing; //Smoothing over frequency
	float* masking; //Activate masking threshold
	float* enable; //For soft bypass (click free bypass)
	// float* attack; //attack time
	// float* release; //release time


	//Parameters values and arrays for the STFT
	int fft_size; //FFTW input size
	int fft_size_2; //FFTW half input size
	int window_combination; // Window combination for the STFT
	float overlap_factor; //oversampling factor for overlap calculations
	float overlap_scale_factor; //Scaling factor for conserving the final amplitude
	int hop; // Hop size for the STFT
	float* window_input; // Input Window values
	float* window_output; // Input Window values
	float* window_count; //Count windows for mean computing
	float max_float; //Auxiliary variable to store FLT_MAX and avoid warning
	float tau; //time constant for soft bypass
	float wet_dry_target; //softbypass target for softbypass
	float wet_dry; //softbypass
	float reduction_coeff; //Gain to apply to the residual noise

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

	//Arrays and variables for getting bins info
	float real_p,imag_n,mag,p2;
	float* fft_magnitude;//magnitude
	float* fft_magnitude_smooth;//magnitude
	float* fft_magnitude_prev;//magnitude spectrum of the previous frame
	float* fft_p2;//power spectrum
	float* fft_p2_smooth;//power spectrum
	float* fft_p2_prev;//power spectum of previous frame
	float* noise_thresholds;
	bool noise_thresholds_availables;

	float* Gk; //gain to be applied
	float* Gk_prev; //past gain applied

	//Loizou algorithm
	float* auto_thresholds; //Reference threshold for louizou algorithm
	float* prev_noise_thresholds;
	float* s_pow_spec;
	float* prev_s_pow_spec;
	float* p_min;
	float* prev_p_min;
	float* speech_p_p;
	float* prev_speech_p_p;

	//masking thresholds
	float* alpha;
	float* beta;
	float* bark_z;
	float max_masked;
	float min_masked;

	// //envelope Follower
	// float* envelope;
	// float attack_coeff;
	// float release_coeff;

	// clock_t start, end;
	// double cpu_time_used;
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
	nrepel->fft_size = FFT_SIZE;
	nrepel->window_combination = WINDOW_COMBINATION;
	nrepel->overlap_factor = OVERLAP_FACTOR;
	nrepel->max_float = FLT_MAX;
	*(nrepel->window_count) = 0.f;
	nrepel->tau = (1.f - exp (-2.f * M_PI * 25.f * 64.f  / nrepel->samp_rate));
	nrepel->noise_thresholds_availables = false;

	nrepel->fft_size_2 = nrepel->fft_size/2;
	nrepel->hop = nrepel->fft_size/nrepel->overlap_factor;
	nrepel->input_latency = nrepel->fft_size - nrepel->hop;
	nrepel->read_ptr = nrepel->input_latency; //the initial position because we are that many samples ahead

	nrepel->in_fifo = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->out_fifo = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_accum = (float*)calloc(nrepel->fft_size,sizeof(float));

	nrepel->window_input = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->window_output = (float*)calloc(nrepel->fft_size,sizeof(float));
	//Window combination computing
	switch(nrepel->window_combination){
		case 0: // HANN-HANN
			fft_window(nrepel->window_input,nrepel->fft_size,0); //STFT input window
			fft_window(nrepel->window_output,nrepel->fft_size,0); //STFT output window
			nrepel->overlap_scale_factor = HANN_HANN_SCALING;
			break;
		case 1: //HAMMING-HANN
			fft_window(nrepel->window_input,nrepel->fft_size,1); //STFT input window
			fft_window(nrepel->window_output,nrepel->fft_size,0); //STFT output window
			nrepel->overlap_scale_factor = HAMMING_HANN_SCALING;
			break;
		case 2: //BLACKMAN-HANN
			fft_window(nrepel->window_input,nrepel->fft_size,2); //STFT input window
			fft_window(nrepel->window_output,nrepel->fft_size,0); //STFT output window
			nrepel->overlap_scale_factor = BLACKMAN_HANN_SCALING;
			break;

	}

	nrepel->input_fft_buffer = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_fft_buffer = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->forward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer, nrepel->output_fft_buffer, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->backward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->output_fft_buffer, nrepel->input_fft_buffer, FFTW_HC2R, FFTW_ESTIMATE);

	nrepel->fft_magnitude = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_magnitude_smooth = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_magnitude_prev = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_p2 = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_p2_smooth = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_p2_prev = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->noise_thresholds = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->auto_thresholds = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	//This was experimentally obteined in louizou paper
	int LF = Freq2Index(1000.f,nrepel->samp_rate,nrepel->fft_size);//1kHz
	int MF = Freq2Index(3000.f,nrepel->samp_rate,nrepel->fft_size);//3kHz
	for (int k = 0;k <= nrepel->fft_size_2; k++){
		if(k < LF){
			nrepel->auto_thresholds[k] = 2.f;
		}
		if(k > LF && k < MF){
			nrepel->auto_thresholds[k] = 2.f;
		}
		if(k > MF){
			nrepel->auto_thresholds[k] = 5.f;
		}
	}

	nrepel->Gk = (float*)malloc((nrepel->fft_size_2+1)*sizeof(float));
	memset(nrepel->Gk, 1, (nrepel->fft_size_2+1)*sizeof(float));
	nrepel->Gk_prev = (float*)malloc((nrepel->fft_size_2+1)*sizeof(float));
	memset(nrepel->Gk_prev, 1, (nrepel->fft_size_2+1)*sizeof(float));

  nrepel->prev_noise_thresholds = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->s_pow_spec = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_s_pow_spec = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->p_min = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_p_min = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->speech_p_p = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_speech_p_p = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	//MASKING THRESHOLDS
	nrepel->alpha = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->beta = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->bark_z = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	memset(nrepel->alpha, 1, (nrepel->fft_size_2+1)*sizeof(float));
	nrepel->max_masked = FLT_MIN;
	nrepel->min_masked = FLT_MAX;
	compute_bark_z(nrepel->bark_z,nrepel->fft_size_2,nrepel->samp_rate);

	// nrepel->envelope = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	return (LV2_Handle)nrepel;
}



static void
connect_port(LV2_Handle instance,
						 uint32_t   port,
						 void*      data) {
	Nrepel* nrepel = (Nrepel*)instance;

	switch ((PortIndex)port) {
		case NREPEL_CAPTURE:
		nrepel->capture_state = (float*)data;
		break;
		case NREPEL_N_AUTO:
		nrepel->auto_state = (float*)data;
		break;
		case NREPEL_AMOUNT:
		nrepel->amount_of_reduction = (float*)data;
		break;
		case NREPEL_STRENGTH:
		nrepel->reduction_strenght = (float*)data;
		break;
		case NREPEL_SMOOTHING:
		nrepel->time_smoothing = (float*)data;
		break;
		case NREPEL_FREQUENCY_SMOOTHING:
		nrepel->frequency_smoothing = (float*)data;
		break;
		// case NREPEL_ATTACK:
		// nrepel->attack = (float*)data;
		// break;
		// case NREPEL_RELEASE:
		// nrepel->release = (float*)data;
		// break;
		case NREPEL_MASKING:
		nrepel->masking = (float*)data;
		break;
		case NREPEL_LATENCY:
		nrepel->report_latency = (float*)data;
		break;
		case NREPEL_WHITENING:
		nrepel->residual_whitening = (float*)data;
		break;
		case NREPEL_RESET:
		nrepel->reset_print = (float*)data;
		break;
		case NREPEL_NOISE_LISTEN:
		nrepel->noise_listen = (float*)data;
		break;
		case NREPEL_ENABLE:
		nrepel->enable = (float*)data;
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

	// //Time execution measurement
	// nrepel->start = clock();
	// //--------------

	//handy variables
	int k;
	unsigned int pos;

	//Inform latency at run call
	*(nrepel->report_latency) = (float) nrepel->input_latency;

	//Softbypass targets in case of disabled or enabled
	if(*(nrepel->enable) == 0.f){ //if disabled
		nrepel->wet_dry_target = 0.f;
	} else { //if enabled
		nrepel->wet_dry_target = 1.f;
	}

	//Interpolate parameters over time softly to bypass without clicks or pops
	nrepel->wet_dry += nrepel->tau * (nrepel->wet_dry_target - nrepel->wet_dry) + FLT_MIN;

	//Reset button state (if on)
	if (*(nrepel->reset_print) == 1.f) {
		memset(nrepel->noise_thresholds, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->Gk, 1, (nrepel->fft_size_2+1)*sizeof(float));
		*(nrepel->window_count) = 0.f;

		memset(nrepel->prev_noise_thresholds, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->s_pow_spec, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_s_pow_spec, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->p_min, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_p_min, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->speech_p_p, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_speech_p_p, 0, (nrepel->fft_size_2+1)*sizeof(float));

		memset(nrepel->alpha, 1, (nrepel->fft_size_2+1)*sizeof(float));
		//memset(nrepel->beta, 0, (nrepel->fft_size_2+1)*sizeof(float));

		nrepel->max_masked = FLT_MIN;
		nrepel->min_masked = FLT_MAX;

		nrepel->noise_thresholds_availables = false;
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

			//Apply windowing
			for (k = 0; k < nrepel->fft_size; k++){
				nrepel->input_fft_buffer[k] = nrepel->in_fifo[k] * nrepel->window_input[k];
			}

			//----------FFT Analysis------------

			//Do transform
			fftwf_execute(nrepel->forward);

			//-----------GET INFO FROM BINS--------------

			//Normalize values to obtain correct magnitude and power values
			for (k = 0; k < nrepel->fft_size; k++){
				nrepel->output_fft_buffer[k] /= nrepel->fft_size;
			}

			//Get the positive spectrum and compute the magnitude
			for (k = 0; k <= nrepel->fft_size_2; k++){
				//Get the half complex spectrum reals and complex
				nrepel->real_p = nrepel->output_fft_buffer[k];
				nrepel->imag_n = nrepel->output_fft_buffer[nrepel->fft_size-k];

				//Get the magnitude and power spectrum
				if(k < nrepel->fft_size){
					nrepel->p2 = (nrepel->real_p*nrepel->real_p + nrepel->imag_n*nrepel->imag_n);
					nrepel->mag = sqrtf(nrepel->p2);//sqrt(real^2+imag^2)
				} else {
					//Nyquist - this is due to half complex transform look at http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html
					nrepel->p2 = nrepel->real_p*nrepel->real_p;
					nrepel->mag = nrepel->real_p;
				}

				//Store values in magnitude and power arrays
				nrepel->fft_p2_prev[k] = nrepel->fft_p2[k]; //store previous value for smoothing
				nrepel->fft_magnitude_prev[k] = nrepel->fft_magnitude[k]; //store previous value for smoothing
				nrepel->fft_p2[k] = sanitize_denormal(nrepel->p2);
				nrepel->fft_magnitude[k] = sanitize_denormal(nrepel->mag);

			}

			/// ---------------PROCESSING--------------------
			//Get noise thresholds if capture is on either auto or manual
			if (*(nrepel->auto_state) == 1.f || *(nrepel->capture_state) == 1.f){
				get_noise_thresholds(*(nrepel->auto_state),
				                     *(nrepel->capture_state),
				                     nrepel->fft_size_2,
				                     nrepel->fft_p2,
				                     nrepel->fft_magnitude,
				                     nrepel->noise_thresholds,
				                     nrepel->auto_thresholds,
				                     nrepel->prev_noise_thresholds,
				                     nrepel->s_pow_spec,
				                     nrepel->prev_s_pow_spec,
				                     nrepel->p_min,
				                     nrepel->prev_p_min,
				                     nrepel->speech_p_p,
				                     nrepel->prev_speech_p_p,
				                     nrepel->window_count);

				nrepel->noise_thresholds_availables = true;
			}

			//Gain calculations and application
			if (nrepel->noise_thresholds_availables){
				//Spectrum preprocessing
				spectral_pre_processing(nrepel->bark_z,
				                    		nrepel->fft_p2,
				                    		nrepel->fft_p2_prev,
				                    		nrepel->fft_p2_smooth,
				                    		nrepel->fft_magnitude,
				                    		nrepel->fft_magnitude_prev,
				                    		nrepel->fft_magnitude_smooth,
				                    		*(nrepel->time_smoothing),
				                    		nrepel->noise_thresholds,
				                    		nrepel->fft_size_2,
				                    		nrepel->alpha,
				                    		//nrepel->beta,
				                    		nrepel->max_masked,
				                    		nrepel->min_masked,
				                    		*(nrepel->reduction_strenght),
				                    		nrepel->Gk,
				                    		nrepel->Gk_prev,
				                    		*(nrepel->masking),
				                    		*(nrepel->frequency_smoothing));

				//Gain application
				gain_application(*(nrepel->amount_of_reduction),
				                 nrepel->fft_size_2,
				                 nrepel->fft_size,
				                 nrepel->output_fft_buffer,
				                 nrepel->Gk,
				                 nrepel->wet_dry,
				                 *(nrepel->residual_whitening),
				                 *(nrepel->noise_listen));

			}
			/// ------------------------------------------------------


			//------------FFT Synthesis-------------

			//Do inverse transform
			fftwf_execute(nrepel->backward);

			//------------OVERLAPADD-------------

			//Accumulate (Overlapadd)
			for(k = 0; k < nrepel->fft_size; k++){
				nrepel->output_accum[k] += nrepel->window_output[k]*nrepel->input_fft_buffer[k]/(float)(nrepel->overlap_scale_factor*nrepel->overlap_factor);
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
