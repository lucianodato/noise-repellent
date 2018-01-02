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
* \file stft.c
* \author Luciano Dato
* \brief Contains a very basic STFT transform abstraction
*/

#include <fftw3.h>
#include "spectral_processing.c"
#include "circular_buffer.c"

/**
* STFT handling struct.
*/
typedef struct
{
  float fft_size;
  float block_size;
  int window_option_input; //Type of input Window for the STFT
  int window_option_output; //Type of output Window for the STFT
  float overlap_factor; //oversampling factor for overlap calculations
  float overlap_scale_factor; //Scaling factor for conserving the final amplitude
  int hop; //Hop size for the STFT
  int input_latency;
  float* input_window;
  float* output_window;
  circular_buffer* input_circular_buffer;
  circular_buffer* output_circular_buffer;
  float* output_accum;
  float* fft_power;
  float* fft_phase;
  float* fft_magnitude;
  float* input_fft_buffer;
  float* output_fft_buffer;
  fftwf_plan forward;
  fftwf_plan backward;
  float fft_size_2;
} STFT_transform;

//-----------STFT private---------------

void
stft_configure(STFT_transform* instance, int block_size, int fft_size,
               int window_option_input, int window_option_output, int overlap_factor)
{
  //Instance configuration
  instance->block_size = block_size;
  instance->fft_size = fft_size;
  instance->fft_size_2 = instance->fft_size/2;
  instance->window_option_input = window_option_input;
  instance->window_option_output = window_option_output;
  instance->overlap_factor = overlap_factor;
  instance->hop = instance->fft_size/instance->overlap_factor;
  instance->input_latency = instance->fft_size - instance->hop;
}

/**
* Wrapper for getting the pre and post processing windows.
* \param input_window array of the input window values
* \param output_window array of the output window values
* \param block_size size of the window arrays
* \param window_option_input input window option
* \param window_option_output output window option
* \param overlap_scale_factor scaling factor for the OLA for configured window options
*/
void
fft_pre_and_post_window(STFT_transform* instance)
{
  //Input window
  switch(instance->window_option_input)
  {
    case 0: // HANN
      fft_window(instance->input_window, instance->block_size, 0); //STFT input window
      break;
    case 1: //HAMMING
      fft_window(instance->input_window, instance->block_size, 1); //STFT input window
      break;
    case 2: //BLACKMAN
      fft_window(instance->input_window, instance->block_size, 2); //STFT input window
      break;
    case 3: //VORBIS
      fft_window(instance->input_window, instance->block_size, 3); //STFT input window
      break;
  }

  //Output window
  switch(instance->window_option_output)
  {
    case 0: // HANN
      fft_window(instance->output_window, instance->block_size, 0); //STFT input window
      break;
    case 1: //HAMMING
      fft_window(instance->output_window, instance->block_size, 1); //STFT input window
      break;
    case 2: //BLACKMAN
      fft_window(instance->output_window, instance->block_size, 2); //STFT input window
      break;
    case 3: //VORBIS
      fft_window(instance->output_window, instance->block_size, 3); //STFT input window
      break;
  }

  //Scaling necessary for perfect reconstruction using Overlapp Add
  instance->overlap_scale_factor = get_window_scale_factor(instance->input_window,
                                                           instance->output_window,
                                                           instance->block_size);
}

/**
* Gets the magnitude and phase spectrum of the complex spectrum. Takimg into account that
* the half complex fft was used half of the spectrum contains the real part and the other
* the imaginary one. Look at http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html for
* more info. DC bin was treated as suggested in http://www.fftw.org/fftw2_doc/fftw_2.html
* \param fft_power the current power spectrum
* \param fft_magnitude the current magnitude spectrum
* \param fft_phase the current phase spectrum
* \param fft_size_2 half of the fft size
* \param fft_size size of the fft
* \param fft_buffer buffer with the complex spectrum of the fft transform
*/
void
get_info_from_bins(float* fft_power, float* fft_magnitude, float* fft_phase,
									 int fft_size_2, int fft_size, float* fft_buffer)
{
	int k;
	float real_p,imag_n,mag,p2,phase;

	//DC bin
	real_p = fft_buffer[0];
	imag_n = 0.f;

	fft_power[0] = real_p*real_p;
	fft_magnitude[0] = real_p;
	fft_phase[0] = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist

	//Get the rest of positive spectrum and compute the magnitude
	for (k = 1; k <= fft_size_2; k++)
	{
		//Get the half complex spectrum reals and complex
		real_p = fft_buffer[k];
		imag_n = fft_buffer[fft_size - k];

		//Get the magnitude, phase and power spectrum
		if(k < fft_size_2)
		{
			p2 = (real_p*real_p + imag_n*imag_n);
			mag = sqrtf(p2);//sqrt(real^2+imag^2)
			phase = atan2f(real_p, imag_n);
		}
		else
		{
			//Nyquist - this is due to half complex transform
			p2 = real_p*real_p;
			mag = real_p;
			phase = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist
		}
		//Store values in magnitude and power arrays (this stores the positive spectrum only)
		fft_power[k] = p2;
		fft_magnitude[k] = mag; //This is not used but part of the STFT transform for generic use
		fft_phase[k] = phase; //This is not used but part of the STFT transform for generic use
	}
}

void
stft_zeropad(STFT_transform* instance)
{
  int k;
  int number_of_zeros = instance->fft_size - instance->block_size;

  //This adds zeros at the end. The right way should be an equal amount at the sides
  for(k = 0; k = number_of_zeros - 1; k++)
  {
    instance->input_fft_buffer[instance->block_size - 1 + k] = 0.f;
  }
}

void
stft_fft_analysis(STFT_transform* instance)
{
  int k;
  //Adding and windowing the frame input values in the center (zero-phasing)
  for (k = 0; k < instance->fft_size; k++)
  {
    instance->input_fft_buffer[k] *= instance->input_window[k];
  }

  //Do transform
  fftwf_execute(instance->forward);
}

void
stft_fft_synthesis(STFT_transform* instance)
{
  int k;
  //Do inverse transform
  fftwf_execute(instance->backward);

  //Normalizing value
  for (k = 0; k < instance->fft_size; k++)
  {
    instance->input_fft_buffer[k] = instance->input_fft_buffer[k]/instance->fft_size;
  }

  //Windowing and scaling
  for(k = 0; k < instance->fft_size; k++)
  {
    instance->input_fft_buffer[k] = (instance->output_window[k]*instance->input_fft_buffer[k])/(instance->overlap_scale_factor * instance->overlap_factor);
  }
}

void
stft_ola(STFT_transform* instance)
{
  int k;
  //Windowing scaling and accumulation
  for(k = 0; k < instance->fft_size; k++)
  {
    instance->output_accum[k] += instance->input_fft_buffer[k];
  }

  //Output samples up to the hop size
  for (k = 0; k < instance->hop; k++)
  {
    cb_write_one(instance->output_circular_buffer, instance->output_accum[k]);
  }

  //shift FFT accumulator the hop size
  memmove(instance->output_accum, instance->output_accum + instance->hop,
          instance->fft_size*sizeof(float));

  //move input FIFO
  for (k = 0; k < instance->input_latency; k++)
  {
    instance->in_fifo[k] = instance->in_fifo[k+instance->hop];
  }
}

//-----------STFT public---------------

void
stft_reset(STFT_transform* instance)
{
  //Reset all arrays
  initialize_array(instance->input_fft_buffer,0.f,instance->fft_size);
  initialize_array(instance->output_fft_buffer,0.f,instance->fft_size);
  initialize_array(instance->input_window,0.f,instance->fft_size);
  initialize_array(instance->output_window,0.f,instance->fft_size);
  cb_reset(instance->input_circular_buffer);
  cb_reset(instance->output_circular_buffer);
  initialize_array(instance->output_accum,0.f,instance->fft_size*2);
  initialize_array(instance->fft_power,0.f,instance->fft_size_2+1);
  initialize_array(instance->fft_magnitude,0.f,instance->fft_size_2+1);
  initialize_array(instance->fft_phase,0.f,instance->fft_size_2+1);
}

void
stft_free(STFT_transform* instance)
{
  fftwf_free(instance->input_fft_buffer);
  fftwf_free(instance->output_fft_buffer);
	fftwf_destroy_plan(instance->forward);
	fftwf_destroy_plan(instance->backward);
	free(instance->input_window);
	free(instance->output_window);
	cb_free(instance->input_circular_buffer);
	cb_free(instance->output_circular_buffer);
	free(instance->output_accum);
	free(instance->fft_power);
	free(instance->fft_magnitude);
	free(instance->fft_phase);
	free(instance);
}


STFT_transform*
stft_init(int block_size, int fft_size,int window_option_input,
          int window_option_output, int overlap_factor)
{
  //Allocate object
  STFT_transform *instance = (STFT_transform*)malloc(sizeof(STFT_transform));

  stft_configure(instance, block_size, fft_size, window_option_input, window_option_output,
                 overlap_factor);

  //Individual array allocation

  //STFT window related
  instance->input_window = (float*)malloc(instance->fft_size * sizeof(float));
  instance->output_window = (float*)malloc(instance->fft_size * sizeof(float));

  //circular buffers init
  instance->input_circular_buffer = cb_init(instance->block_size);
  instance->output_circular_buffer = cb_init(instance->block_size);

  //buffers for OLA
  instance->output_accum = (float*)malloc((instance->fft_size*2) * sizeof(float));

  //FFTW related
  instance->input_fft_buffer = (float*)fftwf_malloc(instance->fft_size * sizeof(float));
  instance->output_fft_buffer = (float*)fftwf_malloc(instance->fft_size * sizeof(float));
  instance->forward = fftwf_plan_r2r_1d(instance->fft_size, instance->input_fft_buffer,
                                       instance->output_fft_buffer, FFTW_R2HC,
                                       FFTW_ESTIMATE);
  instance->backward = fftwf_plan_r2r_1d(instance->fft_size, instance->output_fft_buffer,
                                        instance->input_fft_buffer, FFTW_HC2R,
                                        FFTW_ESTIMATE);

  //Arrays for getting bins info
  instance->fft_power = (float*)malloc((instance->fft_size_2+1) * sizeof(float));
  instance->fft_magnitude = (float*)malloc((instance->fft_size_2+1) * sizeof(float));
  instance->fft_phase = (float*)malloc((instance->fft_size_2+1) * sizeof(float));

  //Initialize all arrays with zeros
  stft_reset(instance);

  //Window combination initialization (pre processing window post processing window)
  fft_pre_and_post_window(instance);

  return instance;
}

void
stft_fill_buffers(STFT_transform* instance, int n_samples, float* input, float* output)
{

}

void
stft_analysis(STFT_transform* instance)
{
  stft_fft_analysis(instance);

  get_info_from_bins(instance->fft_power, instance->fft_magnitude, instance->fft_phase,
                     instance->fft_size_2, instance->fft_size,
                     instance->output_fft_buffer);
}

void
stft_synthesis(STFT_transform* instance)
{
  stft_fft_synthesis(instance);

  stft_ola(instance);
}

void
stft_get_power_spectrum(float* power_spectrum, STFT_transform* instance)
{
  memcpy(power_spectrum, instance->fft_power, sizeof(float)*(instance->fft_size_2+1));
}

void
stft_get_magnitude_spectrum(float* magnitude_spectrum, STFT_transform* instance)
{
  memcpy(magnitude_spectrum, instance->fft_magnitude, sizeof(float)*(instance->fft_size_2+1));
}
