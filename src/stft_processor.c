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
* \file stft_processor.c
* \author Luciano Dato
* \brief Contains an STFT denoiser abstraction
*/

#include "stft_processor.h"
#include "fft_denoiser.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//Window types
#define HANN_WINDOW 0
#define HAMMING_WINDOW 1
#define BLACKMAN_WINDOW 2
#define VORBIS_WINDOW 3

//STFT default values (Hardcoded for now)
#define FFT_SIZE 2048		 //Size of the fft transform
#define INPUT_WINDOW_TYPE 3	 //Input windows for STFT algorithm
#define OUTPUT_WINDOW_TYPE 3 //Output windows for STFT algorithm
#define OVERLAP_FACTOR 4	 //4 is 75% overlap Values bigger than 4 will rescale correctly (if Vorbis windows is not used)

struct STFTProcessor
{
	int fft_size;
	int half_fft_size;
	fftwf_plan forward;
	fftwf_plan backward;
	int window_option_input;	//Type of input Window for the STFT
	int window_option_output;	//Type of output Window for the STFT
	int overlap_factor;			//oversampling factor for overlap calculations
	float overlap_scale_factor; //Scaling factor for conserving the final amplitude
	int hop;					//Hop size for the STFT
	int input_latency;
	int read_position;
	float *input_window;
	float *output_window;
	float *in_fifo;
	float *out_fifo;
	float *output_accum;
	float *input_fft_buffer;
	float *output_fft_buffer;

	//Spectrum information arrays
	float *power_spectrum;
	float *phase_spectrum;
	float *magnitude_spectrum;

	//FFT processor instance
	FFTDenoiser *fft_denoiser;
};

/**
* blackman window values computing.
* \param k bin number
* \param N fft size
*/
static float
blackman(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.42 - 0.5 * cosf(2.f * M_PI * p) + 0.08 * cosf(4.f * M_PI * p);
}

/**
* hanning window values computing.
* \param k bin number
* \param N fft size
*/
static float hanning(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.5 - 0.5 * cosf(2.f * M_PI * p);
}

/**
* hamming window values computing.
* \param k bin number
* \param N fft size
*/
static float hamming(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.54 - 0.46 * cosf(2.f * M_PI * p);
}

/**
* Vorbis window values computing. It satisfies Princen-Bradley criterion so perfect
* reconstruction could be achieved with 50% overlap when used both in Analysis and
* Synthesis
* \param k bin number
* \param N fft size
*/
static float vorbis(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return sinf(M_PI / 2.f * powf(sinf(M_PI * p), 2.f));
}

/**
* Wrapper to compute windows values.
* \param window array for window values
* \param N fft size
* \param window_type type of window
*/
void fft_window(float *window, int N, int window_type)
{
	int k;
	for (k = 0; k < N; k++)
	{
		switch (window_type)
		{
		case BLACKMAN_WINDOW:
			window[k] = blackman(k, N);
			break;
		case HANN_WINDOW:
			window[k] = hanning(k, N);
			break;
		case HAMMING_WINDOW:
			window[k] = hamming(k, N);
			break;
		case VORBIS_WINDOW:
			window[k] = vorbis(k, N);
			break;
		}
	}
}

/**
* Wrapper for getting the pre and post processing windows and adequate scaling factor.
*/
void stft_processor_pre_and_post_window(STFTProcessor *self)
{
	float sum = 0.f;

	//Input window
	switch (self->window_option_input)
	{
	case 0:												   //HANN
		fft_window(self->input_window, self->fft_size, 0); //STFT input window
		break;
	case 1:												   //HAMMING
		fft_window(self->input_window, self->fft_size, 1); //STFT input window
		break;
	case 2:												   //BLACKMAN
		fft_window(self->input_window, self->fft_size, 2); //STFT input window
		break;
	case 3:												   //VORBIS
		fft_window(self->input_window, self->fft_size, 3); //STFT input window
		break;
	}

	//Output window
	switch (self->window_option_output)
	{
	case 0:													//HANN
		fft_window(self->output_window, self->fft_size, 0); //STFT input window
		break;
	case 1:													//HAMMING
		fft_window(self->output_window, self->fft_size, 1); //STFT input window
		break;
	case 2:													//BLACKMAN
		fft_window(self->output_window, self->fft_size, 2); //STFT input window
		break;
	case 3:													//VORBIS
		fft_window(self->output_window, self->fft_size, 3); //STFT input window
		break;
	}

	//Once windows are initialized we can obtain
	//the scaling necessary for perfect reconstruction using Overlapp Add
	for (int i = 0; i < self->fft_size; i++)
		sum += self->input_window[i] * self->output_window[i];

	self->overlap_scale_factor = (sum / (float)(self->fft_size));
}

/**
* Gets the magnitude and phase spectrum of the complex spectrum. Takimg into account that
* the half complex fft was used half of the spectrum contains the real part the other
* the imaginary. Look at http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html for
* more info. DC bin was treated as suggested in http://www.fftw.org/fftw2_doc/fftw_2.html
* \param fft_p2 the current power spectrum
* \param fft_magnitude the current magnitude spectrum
* \param fft_phase the current phase spectrum
* \param fft_size_2 half of the fft size
* \param fft_size size of the fft
* \param fft_buffer buffer with the complex spectrum of the fft transform
*/
void get_info_from_bins(float *fft_p2, float *fft_magnitude, float *fft_phase,
						int fft_size_2, int fft_size, float *fft_buffer)
{
	int k;
	float real_p, imag_n, mag, p2, phase;

	//DC bin
	real_p = fft_buffer[0];
	imag_n = 0.f;

	fft_p2[0] = real_p * real_p;
	fft_magnitude[0] = real_p;
	fft_phase[0] = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist

	//Get the rest of positive spectrum and compute the magnitude
	for (k = 1; k <= fft_size_2; k++)
	{
		//Get the half complex spectrum reals and complex
		real_p = fft_buffer[k];
		imag_n = fft_buffer[fft_size - k];

		//Get the magnitude, phase and power spectrum
		if (k < fft_size_2)
		{
			p2 = (real_p * real_p + imag_n * imag_n);
			mag = sqrtf(p2); //sqrt(real^2+imag^2)
			phase = atan2f(real_p, imag_n);
		}
		else
		{
			//Nyquist - this is due to half complex transform
			p2 = real_p * real_p;
			mag = real_p;
			phase = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist
		}
		//Store values in magnitude and power arrays (this stores the positive spectrum only)
		fft_p2[k] = p2;
		fft_magnitude[k] = mag; //This is not used but part of the STFT transform for generic use
		fft_phase[k] = phase;	//This is not used but part of the STFT transform for generic use
	}
}

/**
* Does the analysis part of the stft for current block.
*/
void stft_processor_analysis(STFTProcessor *self)
{
	int k;

	//Windowing the frame input values in the center (zero-phasing)
	for (k = 0; k < self->fft_size; k++)
	{
		self->input_fft_buffer[k] *= self->input_window[k];
	}

	//Do transform
	fftwf_execute(self->forward);
}

/**
* Does the synthesis part of the stft for current block and then does the OLA method to
* enable the final output.
*/
void stft_processor_synthesis(STFTProcessor *self)
{
	int k;

	//Do inverse transform
	fftwf_execute(self->backward);

	//Normalizing value
	for (k = 0; k < self->fft_size; k++)
	{
		self->input_fft_buffer[k] = self->input_fft_buffer[k] / self->fft_size;
	}

	//Windowing and scaling
	for (k = 0; k < self->fft_size; k++)
	{
		self->input_fft_buffer[k] = (self->output_window[k] * self->input_fft_buffer[k]) / (self->overlap_scale_factor * self->overlap_factor);
	}

	//OVERLAPP-ADD
	//Accumulation
	for (k = 0; k < self->fft_size; k++)
	{
		self->output_accum[k] += self->input_fft_buffer[k];
	}

	//Output samples up to the hop size
	for (k = 0; k < self->hop; k++)
	{
		self->out_fifo[k] = self->output_accum[k];
	}

	//shift FFT accumulator the hop size
	memmove(self->output_accum, self->output_accum + self->hop,
			self->fft_size * sizeof(float));

	//move input FIFO
	for (k = 0; k < self->input_latency; k++)
	{
		self->in_fifo[k] = self->in_fifo[k + self->hop];
	}
}

/**
* Returns the latency needed to be reported to the host.
*/
int stft_processor_get_latency(STFTProcessor *self)
{
	return self->input_latency;
}

/**
* Runs the STFT processing for the given signal by the host.
*/
void stft_processor_run(STFTProcessor *self, int n_samples, const float *input, float *output,
						int enable, int learn_noise, float whitening_factor, float reduction_amount,
						bool residual_listen, float transient_threshold, float masking_ceiling_limit,
						float release, float noise_rescale)
{
	int k;

	for (k = 0; k < n_samples; k++)
	{
		//Read samples given by the host and write samples to the host
		self->in_fifo[self->read_position] = input[k];
		output[k] = self->out_fifo[self->read_position - self->input_latency];
		self->read_position++;

		if (self->read_position >= self->fft_size)
		{
			//Reset read position
			self->read_position = self->input_latency;

			//Fill the fft buffer
			memcpy(self->input_fft_buffer, self->in_fifo, sizeof(float) * self->fft_size);

			//Do analysis
			stft_processor_analysis(self);

			//First get the power, magnitude and phase spectrum
			get_info_from_bins(self->power_spectrum, self->magnitude_spectrum,
							   self->phase_spectrum, self->half_fft_size,
							   self->fft_size, self->output_fft_buffer);

			//Call processing  with the obtained fft transform
			//when stft analysis is applied fft transform values reside in output_fft_buffer
			fft_denoiser_run(self->fft_denoiser, self->power_spectrum, enable, learn_noise, whitening_factor,
							 reduction_amount, residual_listen, transient_threshold, masking_ceiling_limit,
							 release, noise_rescale);

			//Do synthesis
			stft_processor_synthesis(self);
		}
	}
}

/**
* Initializes all dynamics arrays with zeros.
*/
void stft_processor_reset(STFTProcessor *self)
{
	//Reset all arrays
	memset(self->input_fft_buffer, 0.f, self->fft_size);
	memset(self->output_fft_buffer, 0.f, self->fft_size);
	memset(self->input_window, 0.f, self->fft_size);
	memset(self->output_window, 0.f, self->fft_size);
	memset(self->in_fifo, 0.f, self->fft_size);
	memset(self->out_fifo, 0.f, self->fft_size);
	memset(self->output_accum, 0.f, self->fft_size * 2);
	memset(self->power_spectrum, 0.f, self->half_fft_size + 1);
	memset(self->magnitude_spectrum, 0.f, self->half_fft_size + 1);
	memset(self->phase_spectrum, 0.f, self->half_fft_size + 1);
}

/**
* STFT processor initialization and configuration.
*/
STFTProcessor *stft_processor_initialize(int sample_rate)
{
	//Allocate object
	STFTProcessor *self = (STFTProcessor *)malloc(sizeof(STFTProcessor));

	//self configuration
	self->fft_size = FFT_SIZE;
	self->half_fft_size = self->fft_size / 2;
	self->window_option_input = INPUT_WINDOW_TYPE;
	self->window_option_output = OUTPUT_WINDOW_TYPE;
	self->overlap_factor = OVERLAP_FACTOR;
	self->hop = self->fft_size / self->overlap_factor;
	self->input_latency = self->fft_size - self->hop;
	self->read_position = self->input_latency;

	//Individual array allocation

	//STFT window related
	self->input_window = (float *)malloc(self->fft_size * sizeof(float));
	self->output_window = (float *)malloc(self->fft_size * sizeof(float));

	//fifo buffer init
	self->in_fifo = (float *)malloc(self->fft_size * sizeof(float));
	self->out_fifo = (float *)malloc(self->fft_size * sizeof(float));

	//buffer for OLA
	self->output_accum = (float *)malloc((self->fft_size * 2) * sizeof(float));

	//FFTW related
	self->input_fft_buffer = (float *)fftwf_malloc(self->fft_size * sizeof(float));
	self->output_fft_buffer = (float *)fftwf_malloc(self->fft_size * sizeof(float));
	self->forward = fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer,
									  self->output_fft_buffer, FFTW_R2HC,
									  FFTW_ESTIMATE);
	self->backward = fftwf_plan_r2r_1d(self->fft_size, self->output_fft_buffer,
									   self->input_fft_buffer, FFTW_HC2R,
									   FFTW_ESTIMATE);

	//Arrays for getting bins info
	self->power_spectrum = (float *)malloc((self->half_fft_size + 1) * sizeof(float));
	self->magnitude_spectrum = (float *)malloc((self->half_fft_size + 1) * sizeof(float));
	self->phase_spectrum = (float *)malloc((self->half_fft_size + 1) * sizeof(float));

	//Initialize all arrays with zeros
	stft_processor_reset(self);

	//Window combination initialization (pre processing window post processing window)
	stft_processor_pre_and_post_window(self);

	//Spectral processor related
	self->fft_denoiser = fft_denoiser_initialize(self->fft_size, self->fft_size, sample_rate, self->hop);

	return self;
}

/**
* Free allocated memory.
*/
void stft_processor_free(STFTProcessor *self)
{
	fftwf_free(self->input_fft_buffer);
	fftwf_free(self->output_fft_buffer);
	fftwf_destroy_plan(self->forward);
	fftwf_destroy_plan(self->backward);
	free(self->input_window);
	free(self->output_window);
	free(self->in_fifo);
	free(self->out_fifo);
	free(self->output_accum);
	free(self->power_spectrum);
	free(self->magnitude_spectrum);
	free(self->phase_spectrum);
	fft_denoiser_free(self->fft_denoiser);
	free(self);
}

int getHalfSpectralSize(STFTProcessor *self)
{
	return self->half_fft_size;
}

int getSpectralSize(STFTProcessor *self)
{
	return self->fft_size;
}

void setSpectralSize(STFTProcessor *self, int fft_size)
{
	if (!fft_size)
	{
		self->fft_size = FFT_SIZE;
		self->half_fft_size = self->fft_size / 2;
	}

	self->fft_size = fft_size;
	self->half_fft_size = self->fft_size / 2;
}
