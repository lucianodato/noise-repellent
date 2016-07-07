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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "denoise.c"
#include "nestim.c"
#include "extrafunctions.c"

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"

//Noise capture states
#define MANUAL_CAPTURE_OFF_STATE 0
#define MANUAL_CAPTURE_ON_STATE 1
#define AUTO_CAPTURE_STATE 2

//STFT default values
#define DEFAULT_FFT_SIZE 2048 //This should be an even number (Cooley-Turkey)
#define DEFAULT_WINDOW_SIZE 1555 //This should be smaller than FFT size
#define DEFAULT_HOP_SIZE  floor(DEFAULT_WINDOW_SIZE/2) //%50 overlap

///---------------------------------------------------------------------

//LV2 CODE

typedef enum {
	NREPEL_INPUT  = 0,
	NREPEL_OUTPUT = 1,
	NREPEL_CAPTURE = 2,
	NREPEL_AMOUNT = 3,
  NREPEL_WINDOW_TYPE = 4
  NREPEL_LATENCY = 5,

} PortIndex;

typedef struct {
	const float* input;
	float* output;
	float srate;

  //Parameters for the algorithm (user input)
	const int* captstate;
	const float* amountreduc;

  //Parameters for the STFT
	int* samples_needed;
  int* fft_size; //FFTW input size
  int* window_size;
  int* window_type;
  int* hop;
	float* window;
	int* hWS1,hWS2; //first half and second half window size
	int* pFFTsize; //size of positive spectrum, it includes sample 0

  //Temporary buffer for processing
  float* tmpbuf;
	float* current_proc_frame;

  //FFTW related variables
  int* input_size;
  int* output_size;
  float* input_fft_buffer;
  fftwf_complex* output_fft_buffer;
  int flags;
  fftw_plan forward;
  fftw_plan backward;
  float* fft_magnitude;
  float* fft_phase;

	//Store variables
	float* nrepel->noise_print;

} Nrepel;

static LV2_Handle
instantiate(const LV2_Descriptor*     descriptor,
			double                    rate,
			const char*               bundle_path,
			const LV2_Feature* const* features)
{
	Nrepel* nrepel = (Nrepel*)malloc(sizeof(Nrepel));

  //Initialize variables
	nrepel->srate = rate;
  nrepel->fft_size = DEFAULT_FFT_SIZE;
	nrepel->pFFTsize = (nrepel->fft_size/2)+1;
  nrepel->window_size = DEFAULT_WINDOW_SIZE;
	nrepel->window = (float*)fftwf_malloc(sizeof(float)*nrepel->window_size);
	fft_window(nrepel->window,nrepel->window_size,nrepel->window_type);
	nrepel->hop = DEFAULT_HOP_SIZE;
	nrepel->hWS1 = int(floor((nrepel->window_size+1)/2)); //half analysis window size by rounding
	nrepel->hWS2 = int(floor(M/2));//half analysis window size by floor
	nrepel->samples_needed_tmpbfr = nrepel->hWS1 + nrepel->hWS2;
	nrepel->current_proc_frame = (float*)fftwf_malloc(sizeof(float)*nrepel->window_size);

  nrepel->flags = FFTW_ESTIMATE;
  nrepel->input_size = fft_size;
  nrepel->output_size = nrepel->pFFTsize;
  nrepel->input_fft_buffer = (float*)fftwf_malloc(sizeof(float)*fft_size);
  nrepel->output_fft_buffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*output_size);
  nrepel->forward = fftw_plan_dft_r2c_1d(input_size, input_fft_buffer, output_fft_buffer, flags);
  nrepel->backward = fftw_plan_dft_c2r_1d(input_size, output_fft_buffer, input_fft_buffer, flags);

	nrepel->noise_print = (float*)fftwf_malloc(sizeof(float)*nrepel->fft_size);

	return (LV2_Handle)nrepel;
}

static void
connect_port(LV2_Handle instance,
			uint32_t   port,
			void*      data)
{
	Nrepel* nrepel = (Nrepel*)instance;

	switch ((PortIndex)port) {
	case NREPEL_INPUT:
		nrepel->input = (const float*)data;
		break;
	case NREPEL_OUTPUT:
		nrepel->output = (float*)data;
		break;
	case NREPEL_CAPTURE:
		nrepel->captstate = (const int*)data;
		break;
	case NREPEL_AMOUNT:
		nrepel->amountreduc = (const float*)data;
		break;
  case NREPEL_WINDOW_TYPE:
		nrepel->window_type = (const float*)data;
		break;
	case NREPEL_LATENCY:
		nrepel->latency = (int*)data;
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

	const float* input  = nrepel->input;
	float* const output = nrepel->output;

  int type_noise_estimation = nrepel->captstate;

	//-------------------STFT start----------------------
	/*Fill with window_size zeros at the beginning of the received buffer
			and at the end too. This is for correct STFT windowing at the beginning
			and at the end
	*/
	//Reserve the necessary buffer for correct STFT
	nrepel->samples_needed_tmpbfr += n_samples;
  nrepel->tmpbuf = (float*)fftwf_malloc(sizeof(float)*nrepel->samples_needed_tmpbfr);

	// Fill temporary buffer with zeros
	for(int k = 0;k < nrepel->samples_needed_tmpbfr; k++){
		nrepel->tmpbuf[k] = 0;
	}

	//Copy the received signal to the corresponding place in the buffer
	// keep zeros at beginning to center first window at sample 0
	// keep zeros at the end to analyze last sample
	for (int pos = 0; pos < n_samples; pos++) {
		nrepel->tmpbuf[nrepel->hWS2+pos] = input[pos];
	}

	//Auxiliary variables
	int inptr = nrepel->hWS1 //initialize sound pointer in middle of analysis window
	int endptr= nrepel->samples_needed_tmpbfr - nrepel->hWS2 //last sample to start a frame

  //Cycle through the temp buffer given by the host (your daw)
	while (inptr<=endptr) {

		//Get the current frame
		int indx = 0;
		for (int j = inptr-nrepel->hWS1;j<=inptr+nrepel->hWS2;j++){
			nrepel->current_proc_frame[indx] = nrepel->tmpbuf[j];
			indx++;
		}

		//------------------FFTW start-------------------------
		//Do Windowing
		for(int k = 0;k < nrepel->window_size; k++){
			nrepel->current_proc_frame[k] *= nrepel->window[k];
		}

		//initialize the input_fft_buffer with zeros (zeropad)
		for(int k = 0;k < nrepel->input_size; k++){
			nrepel->input_fft_buffer[k] = 0;
		}

		//Zero-Phase Window in fft buffer to avoid phase distortion
		for(int k = 0;k < nrepel->hWS2; k++){
			nrepel->input_fft_buffer[k-1] = nrepel->current_proc_frame[nrepel->output_size-k];
			nrepel->input_fft_buffer[nrepel->hWS1+k] = nrepel->current_proc_frame[k]
		}

		//Do FFT transform
		fftw_execute(forward);

		//Get the magnitude and the phase spectrum
		for (int k = 0; k < nrepel->output_size; k++) {
			//de-interlace FFT buffer
			float real = nrepel->output_fft_buffer[k][0];
			float imag = nrepel->output_fft_buffer[k][1];

			//Magnitude spectrum in dB
			float mag = fabs(real);
			if (mag < FLT_MIN) mag = FLT_MIN;//If mag is smaller than the smallest float assign FLT_MIN
			nrepel->magnitude[k] = 20 * log10f(mag);

			//Phase spectrum
			//Handle denormals
			sanitize_denormal(real);
			sanitize_denormal(imag);
			nrepel->phase[k]=Unwarp(atan2(imag,real),nrepel->output_size);
		}
		//------------------------------------------------------
		//Call denoise function or spectrum estimation function
		switch(nrepel->captstate){
			case MANUAL_CAPTURE_ON_STATE:
				estimate_spectrum(nrepel->magnitude,nrepel->type_noise_estimation,nrepel->noise_print);
				break;
			case MANUAL_CAPTURE_OFF_STATE:
				denoise_signal(nrepel->magnitude);
				break;
			case AUTO_CAPTURE_STATE:
				estimate_spectrum(nrepel->magnitude,nrepel->type_noise_estimation,nrepel->noise_print);
				denoise_signal(nrepel->magnitude);
				break;
		} //switch
		//------------------------------------------------------

		//Reassemble complex spectrum (replace processed magnitude)
		complex aux[(nrepel->output_size-1)*2];
		for (int k = 0; k < nrepel->output_size; k++) {
			// generate positive frequencies
			aux[nrepel->output_size+k] = 10^(nrepel->magnitude[k]/20) * cexp(I*nrepel->phase[k]);
			//generate negative frequencies
			aux[nrepel->output_size-k] = 10^(nrepel->magnitude[nrepel->output_size-k-2]/20) * cexp(-I*nrepel->phase[nrepel->output_size-k-2]);
		}
		//Copy to the output fft buffer to run the plan
		for (int k = 0; k < nrepel->output_size; k++) {
			nrepel->output_fft_buffer[k][0] = crealf(aux[k]);
			nrepel->output_fft_buffer[k][1] = cimagf(aux[k]);
		}

		//Do Inverse FFT
		fftw_execute(backward);

		//------------------FFTW end-------------------------

		//Undo zero-phase window
		for(int k = 0;k < nrepel->hWS2; k++){
			nrepel->current_proc_frame[k] = nrepel->input_fft_buffer[nrepel->hWS1+k];
			nrepel->current_proc_frame[nrepel->output_size-k] = nrepel->input_fft_buffer[k-1]
		}

		//Copy the processed frame to tmpbuf doig overlap-add
		indx = 0;
		for (int j = inptr-nrepel->hWS1;j<=inptr+nrepel->hWS2;j++){
			nrepel->tmpbuf[j] += nrepel->hop * current_proc_frame[indx] = ;
			indx++;
		}
		inptr += nrepel->hop; //advance sound pointer

	}//while

	//-------------------STFT end----------------------

	//Cycle through the processed buffer and output the processed signal
  for (int pos = 0; pos < n_samples; pos++){
    if(nrepel->captstate == MANUAL_CAPTURE_ON_STATE){
      //No processing if the noise spectrum capture state is on
      output[pos] = input[pos];
    }else{
      //Output the processed buffer without added zeros
			output[pos] = nrepel->tmpbuf[pos+nrepel->hWS2];
    }
  }


}

static void
deactivate(LV2_Handle instance)
{
}

static void
cleanup(LV2_Handle instance)
{
	Nrepel* nrepel = (Nrepel*)instance;

	free(nrepel->window);
	free(nrepel->tmpbuf);
	free(nrepel->current_proc_frame)
	free(nrepel->noise_print);
  fftw_free(nrepel->input_fft_buffer);
  fftw_free(nrepel->output_fft_buffer);
  fftw_destroy_plan(nrepel->forward);
  fftw_destroy_plan(nrepel->backward);
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
