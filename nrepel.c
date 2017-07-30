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
#define FFT_SIZE 2048                 //this size should be power of 2
#define WINDOW_TYPE 0          			//0 HANN 1 HAMMING 2 BLACKMAN Input and Output windows for STFT algorithm
#define OVERLAP_FACTOR 4              //4 is 75% overlap Values bigger than 4 will rescale correctly

///---------------------------------------------------------------------

//LV2 CODE

typedef enum {
	NREPEL_AMOUNT = 0,
	NREPEL_NOFFSET = 1,
	NREPEL_RELEASE = 2,
	NREPEL_PF_THRESH = 3,
	NREPEL_WHITENING = 4,
	NREPEL_N_LEARN = 5,
	NREPEL_N_ADAPTIVE = 6,
	NREPEL_RESET = 7,
	NREPEL_NOISE_LISTEN = 8,
	NREPEL_ENABLE = 9,
	NREPEL_LATENCY = 10,
	NREPEL_INPUT = 11,
	NREPEL_OUTPUT = 12,
} PortIndex;

typedef struct {
	const float* input;               //input of samples from host (changing size)
	float* output;                    //output of samples to host (changing size)
	float samp_rate;                  //Sample rate received from the host

	//Parameters for the algorithm (user input)
	float* amount_of_reduction;       //Amount of noise to reduce in dB
	float* noise_thresholds_offset;   //This is to scale the noise profile (over subtraction factor)
	float* pf_threshold;							//threshold for the postfilter detector
	float* release;            	  		//Release time
	float* whitening_factor_pc;				//Whitening amount of the reduction percentage
	float* noise_learn_state;         //Learn Noise state (Manual-Off-Auto)
	float* adaptive_state;            //Autocapture switch
	float* reset_profile;             //Reset Noise switch
	float* noise_listen;              //For noise only listening
	float* enable;                    //For soft bypass (click free bypass)
	float* report_latency;            //Latency necessary

	//Control variables
	bool noise_thresholds_availables; //indicate whether a noise profile is available or no

	//Parameters values and arrays for the STFT
	int fft_size;                     //FFTW input size
	int fft_size_2;                   //FFTW half input size
	int window_option;         				//Type of Window for the STFT
	float overlap_factor;             //oversampling factor for overlap calculations
	float overlap_scale_factor;       //Scaling factor for conserving the final amplitude
	int hop;                          //Hop size for the STFT
	float* window;              			//Window values
	float window_count;               //Count windows for mean computing
	float tau;                        //time constant for soft bypass
	float wet_dry_target;             //softbypass target for softbypass
	float wet_dry;                    //softbypass coeff
	float reduction_coeff;            //Gain to apply to the residual noise
	float release_coeff;							//Release coefficient for Envelopes
	float amount_of_reduction_linear;						//Reduction amount linear value
	float thresholds_offset_linear;		//Threshold offset linear value
	float pf_threshold_linear;				//linear value of the snr threshold
	float whitening_factor;						//Whitening amount of the reduction

	//Buffers for processing and outputting
	int input_latency;
	float* in_fifo;                   //internal input buffer
	float* out_fifo;                  //internal output buffer
	float* output_accum;              //FFT output accumulator
	int read_ptr;                     //buffers read pointer

	//FFTW related arrays
	float* input_fft_buffer;
	float* output_fft_buffer;
	fftwf_plan forward;
	fftwf_plan backward;

	//Postfilter related
	float* input_fft_buffer_ps;
	float* output_fft_buffer_ps;
	fftwf_plan forward_ps;
	fftwf_plan backward_ps;
	float* input_fft_buffer_g;
	float* output_fft_buffer_g;
	fftwf_plan forward_g;
	fftwf_plan backward_g;

	//Arrays and variables for getting bins info
	float real_p,imag_n,mag,p2;
	float* fft_p2;                    //power spectrum
	float* fft_magnitude;             //magnitude spectrum

	//noise related
	float* noise_thresholds_p2;       //captured noise profile power spectrum
	float* noise_thresholds_scaled;   //captured noise profile power spectrum scaled by oversustraction

	//smoothing related
	float* smoothed_spectrum;             //power spectrum to be smoothed
	float* smoothed_spectrum_prev;        //previous frame smoothed power spectrum for envelopes

	//Reduction gains
	float* Gk;			  								//definitive gain

	//Ensemble related
	float* residual_spectrum;
	float* denoised_spectrum;

	//Loizou algorithm
	float* auto_thresholds;           //Reference threshold for louizou algorithm
	float* prev_noise_thresholds;
	float* s_pow_spec;
	float* prev_s_pow_spec;
	float* p_min;
	float* prev_p_min;
	float* speech_p_p;
	float* prev_speech_p_p;

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
instantiate(const LV2_Descriptor*     descriptor,
	    double                    rate,
	    const char*               bundle_path,
	    const LV2_Feature* const* features) {
	//Actual struct declaration
	Nrepel* nrepel = (Nrepel*)malloc(sizeof(Nrepel));

	//Retrieve the URID map callback, and needed URIDs
	for (int i=0; features[i]; ++i) {
		if (!strcmp(features[i]->URI, LV2_URID__map)) {
			nrepel->map = (LV2_URID_Map*)features[i]->data;
		}
	}
	if (!nrepel->map) {
		//bail out: host does not support urid:map
		free(nrepel);
		return NULL;
	}

	nrepel->atom_Vector            = nrepel->map->map(nrepel->map->handle, LV2_ATOM__Vector);
	nrepel->atom_Int               = nrepel->map->map(nrepel->map->handle, LV2_ATOM__Int);
	nrepel->atom_Float             = nrepel->map->map(nrepel->map->handle, LV2_ATOM__Float);
	nrepel->prop_fftsize           = nrepel->map->map(nrepel->map->handle, NREPEL_URI "#fftsize");
	nrepel->prop_nwindow           = nrepel->map->map(nrepel->map->handle, NREPEL_URI "#nwindow");
	nrepel->prop_nrepelFFTp2       = nrepel->map->map(nrepel->map->handle, NREPEL_URI "#FFTp2");

	//Initialize variables
	nrepel->samp_rate = (float)rate;
	nrepel->fft_size = FFT_SIZE;
	nrepel->fft_size_2 = nrepel->fft_size/2;
	nrepel->window_option = WINDOW_TYPE;
	nrepel->overlap_factor = OVERLAP_FACTOR;
	nrepel->hop = nrepel->fft_size/nrepel->overlap_factor;
	nrepel->input_latency = nrepel->fft_size - nrepel->hop;
	nrepel->read_ptr = nrepel->input_latency; //the initial position because we are that many samples ahead
	nrepel->window_count = 0.f;
	nrepel->noise_thresholds_availables = false;
	nrepel->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f  / nrepel->samp_rate));
	nrepel->wet_dry = 0.f;

	nrepel->in_fifo = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->out_fifo = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_accum = (float*)calloc(nrepel->fft_size,sizeof(float));

	nrepel->window = (float*)calloc(nrepel->fft_size,sizeof(float));

	nrepel->input_fft_buffer = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_fft_buffer = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->forward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer, nrepel->output_fft_buffer, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->backward = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->output_fft_buffer, nrepel->input_fft_buffer, FFTW_HC2R, FFTW_ESTIMATE);

	//Postfilter related
	nrepel->input_fft_buffer_ps = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_fft_buffer_ps = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->forward_ps = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer_ps, nrepel->output_fft_buffer_ps, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->backward_ps = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer_ps, nrepel->output_fft_buffer_ps, FFTW_HC2R, FFTW_ESTIMATE);
	nrepel->input_fft_buffer_g = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->output_fft_buffer_g = (float*)calloc(nrepel->fft_size,sizeof(float));
	nrepel->forward_g = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->input_fft_buffer_g, nrepel->output_fft_buffer_g, FFTW_R2HC, FFTW_ESTIMATE);
	nrepel->backward_g = fftwf_plan_r2r_1d(nrepel->fft_size, nrepel->output_fft_buffer_g, nrepel->input_fft_buffer_g, FFTW_HC2R, FFTW_ESTIMATE);

	nrepel->fft_p2 = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->fft_magnitude = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	nrepel->noise_thresholds_p2 = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->noise_thresholds_scaled = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	nrepel->smoothed_spectrum = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->smoothed_spectrum_prev = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	nrepel->Gk = (float*)calloc((nrepel->fft_size),sizeof(float));

	nrepel->residual_spectrum = (float*)calloc((nrepel->fft_size),sizeof(float));
	nrepel->denoised_spectrum = (float*)calloc((nrepel->fft_size),sizeof(float));

	nrepel->auto_thresholds = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_noise_thresholds = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->s_pow_spec = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_s_pow_spec = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->p_min = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_p_min = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->speech_p_p = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));
	nrepel->prev_speech_p_p = (float*)calloc((nrepel->fft_size_2+1),sizeof(float));

	//Window combination initialization (pre processing window post processing window)
	fft_pre_and_post_window(nrepel->window,
				nrepel->fft_size,
				nrepel->window_option,
				&nrepel->overlap_scale_factor);

	//Set initial gain as unity
	memset(nrepel->Gk, 1, (nrepel->fft_size_2+1)*sizeof(float));

	//Compute auto mode initial thresholds
	compute_auto_thresholds(nrepel->auto_thresholds, nrepel->fft_size, nrepel->fft_size_2, nrepel->samp_rate);

	return (LV2_Handle)nrepel;
}



static void
connect_port(LV2_Handle instance,
	     uint32_t   port,
	     void*      data) {
	Nrepel* nrepel = (Nrepel*)instance;

	switch ((PortIndex)port) {
		case NREPEL_AMOUNT:
		nrepel->amount_of_reduction = (float*)data;
		break;
		case NREPEL_NOFFSET:
		nrepel->noise_thresholds_offset = (float*)data;
		break;
		case NREPEL_PF_THRESH:
		nrepel->pf_threshold = (float*)data;
		break;
		case NREPEL_RELEASE:
		nrepel->release = (float*)data;
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
		case NREPEL_NOISE_LISTEN:
		nrepel->noise_listen = (float*)data;
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

	//Parameters values

	/*exponential decay coefficients for envelopes and adaptive noise profiling
		These must take into account the hop size as explained in the following paper
		FFT-BASED DYNAMIC RANGE COMPRESSION*/
	nrepel->release_coeff = expf(-1000.f/(((*(nrepel->release)) * nrepel->samp_rate)/ nrepel->hop) );
	nrepel->amount_of_reduction_linear = from_dB(-1.f * *(nrepel->amount_of_reduction));
	nrepel->thresholds_offset_linear = from_dB(*(nrepel->noise_thresholds_offset));
	nrepel->pf_threshold_linear = from_dB(*(nrepel->pf_threshold));
	nrepel->whitening_factor = *(nrepel->whitening_factor_pc)/100.f;

	//printf("%f\n", nrepel->release_coeff );

	//Reset button state (if on)
	if (*(nrepel->reset_profile) == 1.f) {
		memset(nrepel->noise_thresholds_p2, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->Gk, 1, (nrepel->fft_size_2+1)*sizeof(float));
		nrepel->window_count = 0.f;

		memset(nrepel->prev_noise_thresholds, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->s_pow_spec, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_s_pow_spec, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->p_min, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_p_min, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->speech_p_p, 0, (nrepel->fft_size_2+1)*sizeof(float));
		memset(nrepel->prev_speech_p_p, 0, (nrepel->fft_size_2+1)*sizeof(float));

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
				nrepel->input_fft_buffer[k] = nrepel->in_fifo[k] * nrepel->window[k];
			}

			//----------FFT Analysis------------

			//Do transform
			fftwf_execute(nrepel->forward);

			//-----------GET INFO FROM BINS--------------

			//Get the positive spectrum and compute the magnitude
			for (k = 0; k <= nrepel->fft_size_2; k++){
				//Get the half complex spectrum reals and complex
				nrepel->real_p = nrepel->output_fft_buffer[k];
				nrepel->imag_n = nrepel->output_fft_buffer[nrepel->fft_size-k];

				//Get the magnitude and power spectrum
				if(k < nrepel->fft_size_2){
					nrepel->p2 = (nrepel->real_p*nrepel->real_p + nrepel->imag_n*nrepel->imag_n);
					nrepel->mag = sqrtf(nrepel->p2);//sqrt(real^2+imag^2)
				} else {
					//Nyquist - this is due to half complex transform look at http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html
					nrepel->p2 = nrepel->real_p*nrepel->real_p;
					nrepel->mag = nrepel->real_p;
				}
				//Store values in magnitude and power arrays (this stores the positive spectrum only)
				nrepel->fft_p2[k] = nrepel->p2;
				nrepel->fft_magnitude[k] = nrepel->mag; //This is not used but part of the STFT transform for generic use
			}

			/////////////////////SPECTRAL PROCESSING//////////////////////////

			/*This section countains the specific noise reduction processing blocks
				but it could be replaced with any spectral processing (I'm looking at you future tinkerer)
				Parameters for the STFT transform can be changed at the top of this file
			*/

			//If the spectrum is not silence
			if(!is_empty(nrepel->fft_p2,nrepel->fft_size_2)){
				//If adaptive noise is selected the noise is adapted in time
				if(*(nrepel->adaptive_state) == 1.f) {

					//This has to be revised(issue 8 on github)
					adapt_noise(nrepel->fft_p2,//this is supposed to be the power spectrum in Loizou method
											nrepel->fft_size_2,
											nrepel->noise_thresholds_p2,
											nrepel->auto_thresholds,
											nrepel->prev_noise_thresholds,
											nrepel->s_pow_spec,
											nrepel->prev_s_pow_spec,
											nrepel->p_min,
											nrepel->prev_p_min,
											nrepel->speech_p_p,
											nrepel->prev_speech_p_p);

					nrepel->noise_thresholds_availables = true;
				}

				/*If selected estimate noise spectrum is based on selected portion of signal
				 *do not process the signal
				 */
				if(*(nrepel->noise_learn_state) == 1.f) { //MANUAL
					get_noise_statistics(nrepel->fft_p2,
								nrepel->fft_size_2,
								nrepel->noise_thresholds_p2,
								&nrepel->window_count);

					nrepel->noise_thresholds_availables = true;
				} else {
					//If there is a noise profile reduce noise
					if (nrepel->noise_thresholds_availables == true) {

						//Copy the noise spectrum to the scaled one to be scaled by preprocessing
						memcpy(nrepel->noise_thresholds_scaled,nrepel->noise_thresholds_p2,sizeof(float)*(nrepel->fft_size_2+1));
						//Copy the power spectrum to be smoothed by preprocessing
						memcpy(nrepel->smoothed_spectrum,nrepel->fft_p2,sizeof(float)*(nrepel->fft_size_2+1));

						//Detector smoothing and oversustraction
						preprocessing(nrepel->thresholds_offset_linear,
													nrepel->noise_thresholds_scaled,
													nrepel->smoothed_spectrum,
													nrepel->smoothed_spectrum_prev,
													nrepel->fft_size_2,
													nrepel->Gk,
													nrepel->release_coeff);

						//Supression rule
						spectral_gain(nrepel->smoothed_spectrum,
													nrepel->noise_thresholds_scaled,
													nrepel->fft_size_2,
													*(nrepel->adaptive_state),
													nrepel->Gk);

						//postfilter
						postprocessing(nrepel->fft_size_2,
											    nrepel->fft_size,
													nrepel->fft_p2,
											    nrepel->output_fft_buffer,
													nrepel->input_fft_buffer_ps,
													nrepel->input_fft_buffer_g,
													nrepel->output_fft_buffer_ps,
													nrepel->output_fft_buffer_g,
													&nrepel->forward_g,
													&nrepel->backward_g,
													&nrepel->forward_ps,
													&nrepel->backward_ps,
													nrepel->Gk,
													nrepel->pf_threshold_linear);

						//apply gains
						denoised_calulation(nrepel->fft_size_2,
														    nrepel->fft_size,
														    nrepel->output_fft_buffer,
														    nrepel->denoised_spectrum,
														    nrepel->Gk);

						//residual signal
						residual_calulation(nrepel->fft_size_2,
														    nrepel->fft_size,
														    nrepel->output_fft_buffer,
														    nrepel->residual_spectrum,
														    nrepel->denoised_spectrum,
														    nrepel->Gk,
																nrepel->whitening_factor);

						//Ensemble the final spectrum using residual and denoised
						final_spectrum_ensemble(nrepel->fft_size_2,
														      nrepel->fft_size,
														      nrepel->output_fft_buffer,
														      nrepel->residual_spectrum,
														      nrepel->denoised_spectrum,
																	nrepel->amount_of_reduction_linear,
														      nrepel->wet_dry,
														      *(nrepel->noise_listen));
					}
				}
			}

			///////////////////////////////////////////////////////////

			//------------FFT Synthesis-------------

			//Do inverse transform
			fftwf_execute(nrepel->backward);

			//Windowing and rescaling
			for (k = 0; k < nrepel->fft_size; k++){
				nrepel->input_fft_buffer[k] = (nrepel->window[k]*nrepel->input_fft_buffer[k]) / (nrepel->fft_size * nrepel->overlap_scale_factor * nrepel->overlap_factor);
			}

			//------------OVERLAPADD-------------

			//Accumulate
			for(k = 0; k < nrepel->fft_size; k++){
				nrepel->output_accum[k] += nrepel->input_fft_buffer[k];
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

//spectum struct for noise profile saving
struct FFTVector {
	uint32_t child_size;
	uint32_t child_type;
	float    array[FFT_SIZE/2+1];
};

static LV2_State_Status
savestate(LV2_Handle     instance,
     LV2_State_Store_Function  store,
     LV2_State_Handle          handle,
     uint32_t                  flags,
     const LV2_Feature* const* features)
{
	Nrepel* nrepel = (Nrepel*)instance;

	struct FFTVector vector;

	vector.child_type = nrepel->atom_Float;
	vector.child_size = sizeof(float);

	store(handle, nrepel->prop_fftsize,
			&nrepel->fft_size_2, sizeof(int),
			nrepel->atom_Int, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	store(handle, nrepel->prop_nwindow,
			&nrepel->window_count, sizeof(float),
			nrepel->atom_Float, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	memcpy(vector.array, nrepel->noise_thresholds_p2, sizeof(vector.array));

	store(handle, nrepel->prop_nrepelFFTp2,
			(void*) &vector, sizeof(struct FFTVector),
			nrepel->atom_Vector, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

  return LV2_STATE_SUCCESS;
}

static LV2_State_Status
restorestate(LV2_Handle       instance,
        LV2_State_Retrieve_Function retrieve,
        LV2_State_Handle            handle,
        uint32_t                    flags,
        const LV2_Feature* const*   features)
{
	Nrepel* nrepel = (Nrepel*)instance;
	size_t   size;
	uint32_t type;
	uint32_t valflags;

	//check if state is available
	const int32_t* fftsize = retrieve(handle, nrepel->prop_fftsize, &size, &type, &valflags);
	if (!fftsize || type != nrepel->atom_Int || *fftsize != nrepel->fft_size_2){
		return LV2_STATE_ERR_NO_PROPERTY;
	}

	//check if state is available
	const void* vecFFTp2 = retrieve(handle, nrepel->prop_nrepelFFTp2, &size, &type, &valflags);
	if ( !vecFFTp2 || size != sizeof(struct FFTVector) || type != nrepel->atom_Vector){
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
	if (!strcmp(uri, LV2_STATE__interface)) {
		return &state;
	}
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
