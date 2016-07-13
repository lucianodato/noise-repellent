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


#include "denoise_signal.c"
#include "estimate_spectrum.c"
#include "extra_functions.c"


#include "lv2/lv2plug.in/ns/lv2core/lv2.h"

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"

//Noise capture states
#define MANUAL_CAPTURE_OFF_STATE 0
#define MANUAL_CAPTURE_ON_STATE 1
#define AUTO_CAPTURE_STATE 2

//STFT default values
#define DEFAULT_FFT_SIZE 2048 //This should be an even number (Cooley-Turkey)
#define DEFAULT_WINDOW_TYPE 0 //0 Hann 1 Hamm 2 Black
#define DEFAULT_OVERSAMPLE_FACTOR  2 //2- 50% overlap 4 -75% overlap

///---------------------------------------------------------------------

//LV2 CODE

typedef enum {
	NREPEL_INPUT  = 0,
	NREPEL_OUTPUT = 1,
	NREPEL_CAPTURE = 2,
	NREPEL_AMOUNT = 3,
  NREPEL_LATENCY = 4,
} PortIndex;

typedef struct {
	const float* input;
	float* output;
	float srate; // Sample rate received from the host

  //Parameters for the algorithm (user input)
	int* captstate; // Capture Noise state (Manual-Off-Auto)
	float* amountreduc; // Amount of noise to reduce in dB
	int* report_latency; // Latency necessary

  //Parameters for the STFT
  int fft_size; //FFTW input size
	int fft_size_2; //FFTW half input size
	int windowtype; // Window type for the STFT
	int overlap_factor; //oversampling factor for overlap calculations
	int hop; // Hop size for the STFT
	float* window; // Window values

  //Temporary buffers for processing and outputting
	int input_latency;
  float* inFIFO; //internal input buffer
	float* outFIFO; //internal output buffer
	float* currentFFT; //The frame that is going to be denoised
	float* outputFFTaccum; //FFT output accumulator
	int readPtr; //buffers read pointer

  //FFTW related variables
	//currentFFT is the input of the fftw
	float* input_fft_buffer;
	fftwf_complex* output_fft_buffer;
  int flags;
  fftwf_plan forward;
  fftwf_plan backward;

	float real,imag,mag,phase;
	float* ana_fft_magnitude;//analysis magnitude
  float* ana_fft_phase;//analysis phase
	float* syn_fft_magnitude;//synthesis magnitude
  float* syn_fft_phase;//synthesis phase

	//Store variables
	float* noise_print; // The noise spectrum computed by the captured signal

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
	nrepel->fft_size_2 = (nrepel->fft_size/2);
	nrepel->windowtype = DEFAULT_WINDOW_TYPE;
	nrepel->window = (float*)malloc(sizeof(float)*nrepel->fft_size);
	fft_window(nrepel->window,nrepel->fft_size,nrepel->windowtype); //Fill
	nrepel->overlap_factor = DEFAULT_OVERSAMPLE_FACTOR;
	nrepel->hop = nrepel->fft_size/nrepel->overlap_factor;
	nrepel->input_latency = nrepel->fft_size - nrepel->hop;

	nrepel->inFIFO = (float*)malloc(sizeof(float)*nrepel->fft_size);
	nrepel->outFIFO = (float*)malloc(sizeof(float)*nrepel->fft_size);
	nrepel->outputFFTaccum = (float*)malloc(sizeof(float)*nrepel->fft_size);
	//the initial position because we are that many samples ahead
	nrepel->readPtr = nrepel->input_latency;

  nrepel->flags = FFTW_ESTIMATE;
	nrepel->input_fft_buffer = (float*)malloc(sizeof(float)*nrepel->fft_size);
	nrepel->output_fft_buffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nrepel->fft_size);
	nrepel->forward = fftwf_plan_dft_r2c_1d(nrepel->fft_size, nrepel->input_fft_buffer, nrepel->output_fft_buffer, nrepel->flags);
  nrepel->backward = fftwf_plan_dft_c2r_1d(nrepel->fft_size, nrepel->output_fft_buffer, nrepel->input_fft_buffer, nrepel->flags);

	nrepel->ana_fft_magnitude = (float*)malloc(sizeof(float)*nrepel->fft_size_2);
	nrepel->ana_fft_phase = (float*)malloc(sizeof(float)*nrepel->fft_size_2);
	nrepel->syn_fft_magnitude = (float*)malloc(sizeof(float)*nrepel->fft_size_2);
	nrepel->syn_fft_phase = (float*)malloc(sizeof(float)*nrepel->fft_size_2);

	nrepel->noise_print = (float*)malloc(sizeof(float)*nrepel->fft_size_2);

	//Here we instantiate arrays
	memset(nrepel->inFIFO, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->outFIFO, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->input_fft_buffer, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->outputFFTaccum, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->ana_fft_magnitude, 0, nrepel->fft_size_2*sizeof(float));
	memset(nrepel->ana_fft_phase, 0, nrepel->fft_size_2*sizeof(float));

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
		nrepel->captstate = (int*)data;
		break;
	case NREPEL_AMOUNT:
		nrepel->amountreduc = (float*)data;
		break;
	case NREPEL_LATENCY:
		nrepel->report_latency = (int*)data;
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

	//handy variables
	int k;
	unsigned int pos;

	//Inform latency
	nrepel->report_latency = &nrepel->input_latency;

	//Clear buffers

	//main loop for processing
	for (pos = 0; pos < n_samples; pos++) {
		//Store samples int the input buffer
		nrepel->inFIFO[nrepel->readPtr] = input[pos];
		nrepel->readPtr++;
		//Output samples in the output buffer (even zeros introduced by latency)
		output[pos] = nrepel->inFIFO[nrepel->readPtr-nrepel->input_latency];

		//Once the buffer is full we must process
		if (nrepel->readPtr >= nrepel->fft_size){
			//Reset the input buffer position
			nrepel->readPtr = nrepel->input_latency;

			//Apply windowing
			for (k = 0; k < nrepel->fft_size;k++) {
				nrepel->input_fft_buffer[k] = nrepel->inFIFO[k] * nrepel->window[k];
			}

			//----------FFT Analysis------------
			//Do transform
			fftwf_execute(nrepel->forward);

			//Get the positive spectrum and compute magnitude and phase response
			for (k = 0; k <= nrepel->fft_size_2; k++) {
				nrepel->real = nrepel->output_fft_buffer[k][0];
				nrepel->imag = nrepel->output_fft_buffer[k][1];

				//Calc mag and phase
				nrepel->mag = 2*sqrt(nrepel->real*nrepel->real + nrepel->imag*nrepel->imag);
				nrepel->phase = atan2(nrepel->imag,nrepel->real);

				//Store values in magnitude (in dB) and phase arrays
				nrepel->ana_fft_magnitude[k] = to_dB(nrepel->mag);
				nrepel->ana_fft_phase[k] = nrepel->phase;
			}

			//------------Processing---------------

			//Initialize synthesis arrays and store processing results in them
			memset(nrepel->syn_fft_magnitude, 0, nrepel->fft_size_2*sizeof(float));
			memset(nrepel->syn_fft_phase, 0, nrepel->fft_size_2*sizeof(float));

			//------------FFT Synthesis-------------

			for (k = 0; k <= nrepel->fft_size_2; k++) {
				nrepel->mag = nrepel->syn_fft_magnitude[k];
				nrepel->phase = nrepel->syn_fft_phase[k];

				nrepel->real = from_dB(nrepel->mag)*cos(nrepel->phase);
				nrepel->imag = from_dB(nrepel->mag)*sin(nrepel->phase);

				//Store values in the FFT vector by suming real and the imag part
				nrepel->output_fft_buffer[k][0] = nrepel->real;
				nrepel->output_fft_buffer[k][1] = nrepel->imag;
			}

			//zero negative frequencies of the spectrum
			for (k = nrepel->fft_size_2+1; k < nrepel->fft_size; k++){
				nrepel->output_fft_buffer[k][0] = 0;
			}

			//Do inverse transform
			fftwf_execute(nrepel->backward);

			//Windowing and add to outputFFTaccum
			for(k = 0; k < nrepel->fft_size; k++) {
				nrepel->outputFFTaccum[k] += nrepel->window[k]*nrepel->input_fft_buffer[k]/(nrepel->fft_size_2*nrepel->overlap_factor);
			}
			//Output samples up to the hop size
			for (k = 0; k < nrepel->hop; k++){
				nrepel->outFIFO[k] = nrepel->outputFFTaccum[k];
			}

			//shift FFT accumulator the hop size
			memmove(nrepel->outputFFTaccum, nrepel->outputFFTaccum+nrepel->hop, nrepel->fft_size*sizeof(float));

			//move inputFIFO
			for (k = 0; k < nrepel->input_latency; k++){
				nrepel->inFIFO[k] = nrepel->inFIFO[k+nrepel->hop];
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
	Nrepel* nrepel = (Nrepel*)instance;

	free(nrepel->window);
	free(nrepel->noise_print);
	free(nrepel->inFIFO);
	free(nrepel->outFIFO);
	free(nrepel->outputFFTaccum);
  fftwf_free(nrepel->input_fft_buffer);
  fftwf_free(nrepel->output_fft_buffer);
  fftwf_destroy_plan(nrepel->forward);
  fftwf_destroy_plan(nrepel->backward);
	free(nrepel->ana_fft_magnitude);
	free(nrepel->ana_fft_phase);
	free(nrepel->syn_fft_magnitude);
	free(nrepel->syn_fft_phase);
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
