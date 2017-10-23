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
  float* in_fifo;
  float* out_fifo;
  int read_ptr; //buffers read pointer
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
  instance->read_ptr = instance->input_latency; //the initial position because we are that many samples ahead
}

/**
* Wrapper for getting the pre and post processing windows.
* \param input_window array of the input window values
* \param output_window array of the output window values
* \param frame_size size of the window arrays
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
      fft_window(instance->input_window, instance->frame_size, 0); //STFT input window
      break;
    case 1: //HAMMING
      fft_window(instance->input_window, instance->frame_size, 1); //STFT input window
      break;
    case 2: //BLACKMAN
      fft_window(instance->input_window, instance->frame_size, 2); //STFT input window
      break;
    case 3: //VORBIS
      fft_window(instance->input_window, instance->frame_size, 3); //STFT input window
      break;
  }

  //Output window
  switch(instance->window_option_output)
  {
    case 0: // HANN
      fft_window(instance->output_window, instance->frame_size, 0); //STFT input window
      break;
    case 1: //HAMMING
      fft_window(instance->output_window, instance->frame_size, 1); //STFT input window
      break;
    case 2: //BLACKMAN
      fft_window(instance->output_window, instance->frame_size, 2); //STFT input window
      break;
    case 3: //VORBIS
      fft_window(instance->output_window, instance->frame_size, 3); //STFT input window
      break;
  }

  //Scaling necessary for perfect reconstruction using Overlapp Add
  instance->overlap_scale_factor = get_window_scale_factor(instance->input_window,
                                                           instance->output_window,
                                                           instance->frame_size);
}

/**
* Gets the magnitude and phase spectrum of the complex spectrum. Takimg into account that
* the half complex fft was used half of the spectrum contains the real part the other
* the imaginary. Look at http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html for
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

  //----------FFT Analysis------------

  //Do transform
  fftwf_execute(instance->forward);
}

void
stft_fft_synthesis(STFT_transform* instance)
{
  //------------FFT Synthesis-------------

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
  //Windowing scaling and accumulation
  for(k = 0; k < instance->fft_size; k++)
  {
    instance->output_accum[k] += instance->input_fft_buffer[k];
  }

  //Output samples up to the hop size
  for (k = 0; k < instance->hop; k++)
  {
    instance->out_fifo[k] = instance->output_accum[k];
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
  initialize_array(instance->input_fft_buffer,0.f,self->fft_size);
  initialize_array(instance->output_fft_buffer,0.f,self->fft_size);
  initialize_array(instance->input_window,0.f,self->fft_size);
  initialize_array(instance->output_window,0.f,self->fft_size);
  initialize_array(instance->in_fifo,0.f,self->fft_size);
  initialize_array(instance->out_fifo,0.f,self->fft_size);
  initialize_array(instance->output_accum,0.f,self->fft_size*2);
  initialize_array(instance->fft_power,0.f,self->fft_size_2+1);
  initialize_array(instance->fft_magnitude,0.f,self->fft_size_2+1);
  initialize_array(instance->fft_phase,0.f,self->fft_size_2+1);
}

void
stft_init(STFT_transform* instance, int block_size, int fft_size,int window_option_input,
          int window_option_output, int overlap_factor)
{
  stft_configure(instance, block_size, fft_size, window_option_input, window_option_output,
                 overlap_factor);

  //STFT window related
  instance->input_window = (float*)malloc(instance->fft_size, sizeof(float));
  instance->output_window = (float*)malloc(instance->fft_size, sizeof(float));

  //buffers for OLA
  instance->in_fifo = (float*)malloc(instance->fft_size, sizeof(float));
  instance->out_fifo = (float*)malloc(instance->fft_size, sizeof(float));
  instance->output_accum = (float*)malloc(instance->fft_size*2, sizeof(float));

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
  instance->fft_power = (float*)malloc((instance->fft_size_2+1), sizeof(float));
  instance->fft_magnitude = (float*)malloc((instance->fft_size_2+1), sizeof(float));
  instance->fft_phase = (float*)malloc((instance->fft_size_2+1), sizeof(float));

  //Initialize all arrays with zeros
  stft_reset(instance);

  //Window combination initialization (pre processing window post processing window)
  fft_pre_and_post_window(instance);
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
	free(instance->in_fifo);
	free(instance->out_fifo);
	free(instance->output_accum);
	free(instance->fft_power);
	free(instance->fft_magnitude);
	free(instance->fft_phase);
	free(instance);
}

void
stft_run(STFT_transform* instance, uint32_t n_samples, float* input, float* output)
{
  stft_analysis(instance);

  get_info_from_bins(instance->fft_power, instance->fft_magnitude, instance->fft_phase,
                     instance->fft_size_2, instance->fft_size,
                     instance->output_fft_buffer);

  stft_synthesis(instance);
}

void
stft_get_p2_spectrum(STFT_transform* instance)
{

}
