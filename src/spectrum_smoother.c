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
* \file spectrum_smoother.c
* \author Luciano Dato
* \brief Contains a spectrum smoother abstraction
*/

#include "spectrum_smoother.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct SpectralSmoother
{
	//General parameters
	int fft_size;
	int half_fft_size;
	int samp_rate;
	int hop;

	//Ensemble related
	//Spectrum
	float *noise_spectrum;
	float *signal_spectrum;

	//smoothing related
	float *smoothed_spectrum;	   //power spectrum to be smoothed
	float *smoothed_spectrum_prev; //previous frame smoothed power spectrum for envelopes

	float release_coefficient; //reference smoothing value
};

/*
* Exponential decay coefficients for envelopes and adaptive noise profiling
* These must take into account the hop size as explained in the following paper
* FFT-BASED DYNAMIC RANGE COMPRESSION
*/
void get_release_coefficient(SpectralSmoother *self, float release)
{
	if (release != 0.f) //This allows to turn off smoothing with 0 ms in order to use masking only
	{
		self->release_coefficient = expf(-1000.f / (((release)*self->samp_rate) / self->hop));
	}
	else
	{
		self->release_coefficient = 0.f; //This avoids incorrect results when moving sliders rapidly
	}
}

/**
* Spectral time smoothing by applying a release envelope. This seems to work better than * using time smoothing directly or McAulay & Malpass modification.
*/
void apply_time_envelope(SpectralSmoother *self)
{
	int k;

	for (k = 0; k <= self->half_fft_size; k++)
	{
		//It doesn't make much sense to have an attack slider when there is time smoothing
		if (self->smoothed_spectrum[k] > self->smoothed_spectrum_prev[k])
		{
			//Release (when signal is incrementing in amplitude)
			self->smoothed_spectrum[k] = self->release_coefficient * self->smoothed_spectrum_prev[k] + (1.f - self->release_coefficient) * self->smoothed_spectrum[k];
		}
	}
}

void spectral_smoothing_run(SpectralSmoother *self, float release)
{
	get_release_coefficient(self, release);

	memcpy(self->smoothed_spectrum, self->signal_spectrum, sizeof(float) * (self->half_fft_size + 1));

	apply_time_envelope(self);

	memcpy(self->smoothed_spectrum_prev, self->smoothed_spectrum, sizeof(float) * (self->half_fft_size + 1));
}

/**
* Reset dynamic arrays to zero.
*/
void spectral_smoothing_reset(SpectralSmoother *self)
{
	//Reset all arrays
	memset(self->signal_spectrum, 0.f, self->half_fft_size + 1);
	memset(self->noise_spectrum, 0.f, self->half_fft_size + 1);
	memset(self->smoothed_spectrum, 0.f, self->half_fft_size + 1);
	memset(self->smoothed_spectrum_prev, 0.f, self->half_fft_size + 1);

	self->release_coefficient = 0.f;
}

/**
* Gain estimator initialization and configuration.
*/
SpectralSmoother *spectral_smoothing_initialize(int fft_size, int samp_rate, int hop)
{
	//Allocate object
	SpectralSmoother *self = (SpectralSmoother *)malloc(sizeof(SpectralSmoother));

	//Configuration
	self->fft_size = fft_size;
	self->half_fft_size = self->fft_size / 2;
	self->samp_rate = samp_rate;
	self->hop = hop;

	//spectrum allocation
	self->signal_spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->noise_spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->smoothed_spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->smoothed_spectrum_prev = (float *)calloc((self->half_fft_size + 1), sizeof(float));

	//Reset all values
	spectral_smoothing_reset(self);

	return self;
}

/**
* Free allocated memory.
*/
void spectral_smoothing_free(SpectralSmoother *self)
{
	free(self->noise_spectrum);
	free(self->signal_spectrum);
	free(self->smoothed_spectrum);
	free(self->smoothed_spectrum_prev);
	free(self);
}
