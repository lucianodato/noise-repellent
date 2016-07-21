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
#define DEFAULT_WINDOW_SIZE 2001 //This should be an odd number fft/2<=WS<fft
#define DEFAULT_WINDOW_TYPE 0 //0 Hann 1 Hamming 2 Blackman
#define DEFAULT_OVERLAP_FACTOR 4 //2 is 50% and 4 is 75% overlap

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
	const float* input; //input of samples from host (changing size)
	float* output; //output of samples to host (changing size)
	float samp_rate; // Sample rate received from the host

  //Parameters for the algorithm (user input)
	int* capt_state; // Capture Noise state (Manual-Off-Auto)
	int* fft_option; //selected fft size by the user
	float* amount_reduc; // Amount of noise to reduce in dB
	int* report_latency; // Latency necessary

  //Parameters for the STFT
  int fft_size; //FFTW input size
	int fft_size_2; //FFTW half input size
	int window_size; //Window size odd number
	int w_fh, w_sh; //Window firts half and second half
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
	fftwf_complex* output_fft_buffer;
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
	//Actual struct declaration
	Nrepel* nrepel = (Nrepel*)malloc(sizeof(Nrepel));

  //Initialize variables
	nrepel->samp_rate = rate;
  nrepel->fft_size = DEFAULT_FFT_SIZE;
	nrepel->window_size = DEFAULT_WINDOW_SIZE;
	nrepel->window_type = DEFAULT_WINDOW_TYPE;
	nrepel->overlap_factor = DEFAULT_OVERLAP_FACTOR;

	nrepel->fft_size_2 = nrepel->fft_size/2;
	nrepel->w_fh = floor(nrepel->window_size+1)/2;
	nrepel->w_sh = floor(nrepel->window_size)/2;
	nrepel->hop = nrepel->window_size/nrepel->overlap_factor;
	nrepel->input_latency = nrepel->window_size - nrepel->hop;
	nrepel->read_ptr = nrepel->input_latency; //the initial position because we are that many samples ahead

	nrepel->in_fifo = (float*)malloc(sizeof(float)*nrepel->window_size);
	nrepel->out_fifo = (float*)malloc(sizeof(float)*nrepel->window_size);
	nrepel->output_accum = (float*)malloc(sizeof(float)*nrepel->window_size);

	nrepel->window = (float*)malloc(sizeof(float)*nrepel->window_size);
	nrepel->input_fft_buffer = (float*)malloc(sizeof(float)*nrepel->fft_size);
	nrepel->output_fft_buffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nrepel->fft_size);
	nrepel->forward = fftwf_plan_dft_r2c_1d(nrepel->fft_size, nrepel->input_fft_buffer, nrepel->output_fft_buffer, FFTW_ESTIMATE);
  nrepel->backward = fftwf_plan_dft_c2r_1d(nrepel->fft_size, nrepel->output_fft_buffer, nrepel->input_fft_buffer, FFTW_ESTIMATE);

	nrepel->ana_fft_magnitude = (float*)malloc(sizeof(float)*nrepel->fft_size);
	nrepel->ana_fft_phase = (float*)malloc(sizeof(float)*nrepel->fft_size);
	nrepel->syn_fft_magnitude = (float*)malloc(sizeof(float)*nrepel->fft_size);
	nrepel->syn_fft_phase = (float*)malloc(sizeof(float)*nrepel->fft_size);

	nrepel->noise_print = (float*)malloc(sizeof(float)*nrepel->fft_size);

	//Here we initialize arrays with zeros
	memset(nrepel->in_fifo, 0, nrepel->window_size*sizeof(float));
	memset(nrepel->out_fifo, 0, nrepel->window_size*sizeof(float));
	memset(nrepel->input_fft_buffer, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->output_fft_buffer, 0, nrepel->fft_size*sizeof(fftwf_complex));
	memset(nrepel->output_accum, 0, nrepel->window_size*sizeof(float));
	memset(nrepel->ana_fft_magnitude, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->ana_fft_phase, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->syn_fft_magnitude, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->syn_fft_phase, 0, nrepel->fft_size*sizeof(float));
	memset(nrepel->noise_print, 0, nrepel->fft_size*sizeof(float));

	fft_window(nrepel->window,nrepel->window_size,nrepel->window_type); //Init window

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
		nrepel->capt_state = (int*)data;
		break;
	case NREPEL_AMOUNT:
		nrepel->amount_reduc = (float*)data;
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

	//handy variables
	int k;
	unsigned int pos;

	//Inform latency at run call
	*(nrepel->report_latency) = nrepel->input_latency;

	//main loop for processing
	for (pos = 0; pos < n_samples; pos++){
		//Store samples int the input buffer
		nrepel->in_fifo[nrepel->read_ptr] = nrepel->input[pos];
		//Output samples in the output buffer (even zeros introduced by latency)
		nrepel->output[pos] = nrepel->out_fifo[nrepel->read_ptr - nrepel->input_latency];
		//Now move the read pointer
		nrepel->read_ptr++;

		//Once the buffer is full we can do stuff
		if (nrepel->read_ptr >= nrepel->window_size){
			//Reset the input buffer position
			nrepel->read_ptr = nrepel->input_latency;

			//Reinitialize fft buffer (zeros are already there)
			memset(nrepel->input_fft_buffer, 0, nrepel->fft_size*sizeof(float));

			//Apply windowing (Could be any window type)
			for (k = 0; k < nrepel->window_size; k++){
				nrepel->input_fft_buffer[k] = nrepel->in_fifo[k] * nrepel->window[k];
			}

			//Zerophase window shifting (zeropad the center to avoid phase distortion)
			zero_phase_window(nrepel->input_fft_buffer,nrepel->w_sh,nrepel->window_size);

			//----------FFT Analysis------------
			//Do transform
			fftwf_execute(nrepel->forward);

			//Get the positive spectrum and compute magnitude and phase response
			for (k = nrepel->fft_size_2-1; k <= nrepel->fft_size_2; k++){
				nrepel->real = nrepel->output_fft_buffer[k][0];
				nrepel->imag = nrepel->output_fft_buffer[k][1];

				//Get mag and phase
				nrepel->mag = sanitize_denormal(sqrtf(nrepel->real*nrepel->real + nrepel->imag*nrepel->imag));
				nrepel->phase = sanitize_denormal(atan2f(nrepel->imag,nrepel->real));

				//Store values in magnitude and phase arrays
				nrepel->ana_fft_magnitude[k] = nrepel->mag;
				nrepel->ana_fft_phase[k] = nrepel->phase;
			}

			//------------Processing---------------

			//Call processing functions and send nrepel->ana_fft_magnitude[k]
			for (k = nrepel->fft_size_2-1; k <= nrepel->fft_size_2; k++){
					nrepel->syn_fft_magnitude[k] = nrepel->ana_fft_magnitude[k];
					nrepel->syn_fft_phase[k] = nrepel->ana_fft_phase[k];
			}

			//------------FFT Synthesis-------------

			for (k = nrepel->fft_size_2-1; k <= nrepel->fft_size_2; k++){
				nrepel->mag = nrepel->syn_fft_magnitude[k];
				nrepel->phase = nrepel->syn_fft_phase[k];

				nrepel->real = sanitize_denormal(nrepel->mag*cosf(nrepel->phase));
				nrepel->imag = sanitize_denormal(nrepel->mag*sinf(nrepel->phase));

				//Store values in the FFT vector by suming real and the imag part
				nrepel->output_fft_buffer[k][0] = nrepel->real;
				nrepel->output_fft_buffer[k][1] = nrepel->imag;
			}

			//Make the negative spectrum zero
			for (k = 0; k < nrepel->fft_size_2; k++){
				nrepel->output_fft_buffer[k][0] = 0.f;
				nrepel->output_fft_buffer[k][1] = 0.f;
			}

			//Do inverse transform
			fftwf_execute(nrepel->backward);

			//Undo zerophase shift
			zero_phase_window(nrepel->input_fft_buffer,nrepel->w_fh,nrepel->window_size);

			//Scaling and add to output_accum
			for(k = 0; k < nrepel->window_size; k++){
				//Scaling FFT (because is unnormalized)
				nrepel->input_fft_buffer[k] /= nrepel->fft_size;
				//Accumulate (Overlapadd)
				nrepel->output_accum[k] += nrepel->input_fft_buffer[k]*nrepel->hop;
			}

			//Output samples up to the hop size
			for (k = 0; k < nrepel->hop; k++){
				nrepel->out_fifo[k] = nrepel->output_accum[k];
			}

			//shift FFT accumulator the hop size
			memmove(nrepel->output_accum, nrepel->output_accum + nrepel->hop, nrepel->window_size*sizeof(float));

			for (k = 0; k < nrepel->input_latency; k++){
				nrepel->in_fifo[k] = nrepel->in_fifo[k+nrepel->hop];
			}

			//move inputFIFO
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
	Nrepel* nrepel = (Nrepel*)instance;

	free(nrepel->window);
	free(nrepel->noise_print);
	free(nrepel->in_fifo);
	free(nrepel->out_fifo);
	free(nrepel->output_accum);
	free(nrepel->input_fft_buffer);
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
