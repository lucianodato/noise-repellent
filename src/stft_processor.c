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

#include "stft_processor.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define INPUT_WINDOW_TYPE 3
#define OUTPUT_WINDOW_TYPE 3

typedef enum
{
	HANN_WINDOW = 0,
	HAMMING_WINDOW = 1,
	BLACKMAN_WINDOW = 2,
	VORBIS_WINDOW = 3
} WindowTypes;

struct STFTProcessor
{
	int fft_size;
	int half_fft_size;
	fftwf_plan forward;
	fftwf_plan backward;
	int window_option_input;
	int window_option_output;
	int overlap_factor;
	float overlap_scale_factor;
	int hop;
	int input_latency;
	int read_position;
	float *input_window;
	float *output_window;
	float *in_fifo;
	float *out_fifo;
	float *output_accum;
	float *input_fft_buffer;
	float *output_fft_buffer;

	FFTDenoiser *fft_denoiser;
};

static float blackman(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.42 - 0.5 * cosf(2.f * M_PI * p) + 0.08 * cosf(4.f * M_PI * p);
}

static float hanning(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.5 - 0.5 * cosf(2.f * M_PI * p);
}

static float hamming(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.54 - 0.46 * cosf(2.f * M_PI * p);
}

static float vorbis(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return sinf(M_PI / 2.f * powf(sinf(M_PI * p), 2.f));
}

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

void stft_processor_pre_and_post_window(STFTProcessor *self)
{
	float sum = 0.f;

	//Input window
	switch ((WindowTypes)self->window_option_input)
	{
	case HANN_WINDOW:
		fft_window(self->input_window, self->fft_size, HANN_WINDOW);
		break;
	case HAMMING_WINDOW:
		fft_window(self->input_window, self->fft_size, HAMMING_WINDOW);
		break;
	case BLACKMAN_WINDOW:
		fft_window(self->input_window, self->fft_size, BLACKMAN_WINDOW);
		break;
	case VORBIS_WINDOW:
		fft_window(self->input_window, self->fft_size, VORBIS_WINDOW);
		break;
	}

	//Output window
	switch ((WindowTypes)self->window_option_output)
	{
	case HANN_WINDOW:
		fft_window(self->output_window, self->fft_size, HANN_WINDOW);
		break;
	case HAMMING_WINDOW:
		fft_window(self->output_window, self->fft_size, HAMMING_WINDOW);
		break;
	case BLACKMAN_WINDOW:
		fft_window(self->output_window, self->fft_size, BLACKMAN_WINDOW);
		break;
	case VORBIS_WINDOW:
		fft_window(self->output_window, self->fft_size, VORBIS_WINDOW);
		break;
	}

	for (int i = 0; i < self->fft_size; i++)
		sum += self->input_window[i] * self->output_window[i];

	self->overlap_scale_factor = (sum / (float)(self->fft_size));
}

void stft_processor_analysis(STFTProcessor *self)
{
	int k;

	for (k = 0; k < self->fft_size; k++)
	{
		self->input_fft_buffer[k] *= self->input_window[k];
	}

	fftwf_execute(self->forward);
}

void stft_processor_synthesis(STFTProcessor *self)
{
	int k;

	fftwf_execute(self->backward);

	for (k = 0; k < self->fft_size; k++)
	{
		self->input_fft_buffer[k] = self->input_fft_buffer[k] / self->fft_size;
	}

	for (k = 0; k < self->fft_size; k++)
	{
		self->input_fft_buffer[k] = (self->output_window[k] * self->input_fft_buffer[k]) / (self->overlap_scale_factor * self->overlap_factor);
	}

	for (k = 0; k < self->fft_size; k++)
	{
		self->output_accum[k] += self->input_fft_buffer[k];
	}

	for (k = 0; k < self->hop; k++)
	{
		self->out_fifo[k] = self->output_accum[k];
	}

	memmove(self->output_accum, self->output_accum + self->hop,
			self->fft_size * sizeof(float));

	for (k = 0; k < self->input_latency; k++)
	{
		self->in_fifo[k] = self->in_fifo[k + self->hop];
	}
}

int stft_processor_get_latency(STFTProcessor *self)
{
	return self->input_latency;
}

void stft_processor_run(STFTProcessor *self, NoiseProfile *noise_profile, int n_samples, const float *input, float *output)
{
	for (int k = 0; k < n_samples; k++)
	{
		self->in_fifo[self->read_position] = input[k];
		output[k] = self->out_fifo[self->read_position - self->input_latency];
		self->read_position++;

		if (self->read_position >= self->fft_size)
		{
			self->read_position = self->input_latency;

			memcpy(self->input_fft_buffer, self->in_fifo, sizeof(float) * self->fft_size);

			stft_processor_analysis(self);

			fft_denoiser_run(self->fft_denoiser, noise_profile, self->output_fft_buffer);

			stft_processor_synthesis(self);
		}
	}
}

STFTProcessor *stft_processor_initialize(FFTDenoiser *fft_denoiser, int fft_size, int overlap_factor)
{
	STFTProcessor *self = (STFTProcessor *)calloc(1, sizeof(STFTProcessor));

	set_spectral_size(self, fft_size);
	self->window_option_input = INPUT_WINDOW_TYPE;
	self->window_option_output = OUTPUT_WINDOW_TYPE;
	self->overlap_factor = overlap_factor;
	self->hop = self->fft_size / self->overlap_factor;
	self->input_latency = self->fft_size - self->hop;
	self->read_position = self->input_latency;

	self->input_window = (float *)calloc(self->fft_size, sizeof(float));
	self->output_window = (float *)calloc(self->fft_size, sizeof(float));

	self->in_fifo = (float *)calloc(self->fft_size, sizeof(float));
	self->out_fifo = (float *)calloc(self->fft_size, sizeof(float));

	self->output_accum = (float *)calloc((self->fft_size * 2), sizeof(float));

	self->input_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
	self->output_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
	self->forward = fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer,
									  self->output_fft_buffer, FFTW_R2HC,
									  FFTW_ESTIMATE);
	self->backward = fftwf_plan_r2r_1d(self->fft_size, self->output_fft_buffer,
									   self->input_fft_buffer, FFTW_HC2R,
									   FFTW_ESTIMATE);

	stft_processor_pre_and_post_window(self);

	self->fft_denoiser = fft_denoiser;

	return self;
}

void stft_processor_free(STFTProcessor *self)
{
	free(self->input_fft_buffer);
	free(self->output_fft_buffer);
	fftwf_destroy_plan(self->forward);
	fftwf_destroy_plan(self->backward);
	free(self->input_window);
	free(self->output_window);
	free(self->in_fifo);
	free(self->out_fifo);
	free(self->output_accum);
	free(self);
}

void set_spectral_size(STFTProcessor *self, int fft_size)
{
	self->fft_size = fft_size;
	self->half_fft_size = self->fft_size / 2;
}
