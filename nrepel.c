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
#define DEFAULT_WINDOW_SIZE 1555 //This should be an odd number (zerophase window)

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
  int* fftsize; //FFTW input size
  int* window_size;
  int* window_type;
  int* overlap;

  //Temporary buffer to reach fftsize
  float* tmpbuf;
  int* bufptr; //buffer position pointer

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
  nrepel->bufptr = 0;
  nrepel->tmpbuf = (float*)fftwf_malloc(sizeof(float)*fftsize);
  nrepel->fftsize = DEFAULT_FFT_SIZE;
  nrepel->window_size = DEFAULT_WINDOW_SIZE;

  nrepel->flags = FFTW_ESTIMATE;
  nrepel->input_size = fftsize;
  nrepel->output_size = (nrepel->input_size/2 - 1);
  nrepel->input_fft_buffer = (float*)fftwf_malloc(sizeof(float)*fftsize);
  nrepel->output_fft_buffer = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*output_size);
  nrepel->forward = fftw_plan_dft_r2c_1d(input_size, input_fft_buffer, output_fft_buffer, flags);
  nrepel->backward = fftw_plan_dft_c2r_1d(input_size, output_fft_buffer, input_fft_buffer, flags);

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

  //1-Cycle through the buffer given by the host (your daw)
	for (uint32_t pos = 0; pos < n_samples; pos++) {

      //2-Store those samples in a temporary buffer
      nrepel->tmpbuf[nrepel->bufptr] = input[pos];
      nrepel->bufptr++;

      //3-If the window size selected is reached do processing
      if (nrepel->bufptr > nrepel->window_size) {
				nrepel->bufptr = 0; //restarting the pointer

        //4-Do all STFT stuff before processing

        //5-Zeropad the array if necessary

        //6-Do FFT transform
        fftw_execute(forward);

        //7-Get the magnitude and the phase spectrum
    		for (int k = 0; k <= nrepel->output_size; k++) {
    			//de-interlace FFT buffer
    			float real = nrepel->output_fft_buffer[k][0];
    			float imag = nrepel->output_fft_buffer[k][1];

    			//compute magnitude and phase
    			nrepel->magnitude[k] = 2.*sqrt(real*real + imag*imag);
    			nrepel->phase[k]=atan2(imag,real);
    		}

        //8-Call denoise function or spectrum estimation function
        switch(nrepel->captstate){
          case MANUAL_CAPTURE_ON_STATE:
            estimate_spectrum(nrepel->magnitude,nrepel->type_noise_estimation);
            break;
          case MANUAL_CAPTURE_OFF_STATE:
            denoise_signal(nrepel->magnitude);
            break;
          case AUTO_CAPTURE_STATE:
            estimate_spectrum(nrepel->magnitude,nrepel->type_noise_estimation);
            denoise_signal(nrepel->magnitude);
            break;
        }

        //9-Reassemble complex spectrum (replace processed magnitude)
        for (int k = 0; k <= nrepel->output_size; k++) {
          float real = nrepel->magnitude[k] * cos(nrepel->phase[k]);
    			float imag = nrepel->magnitude[k] * cos(nrepel->phase[k]);

    			nrepel->output_fft_buffer[k][0] = real;
    			nrepel->output_fft_buffer[k][1] = imag;
    		}

        //10-Do Inverse FFT
        fftw_execute(backward);

        //11-Assign the processed samples to the output buffer

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
  fftw_free(input_fft_buffer);
  fftw_free(output_fft_buffer);
  fftw_destroy_plan(forward);
  fftw_destroy_plan(backward);
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
