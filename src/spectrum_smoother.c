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

#include "extra_functions.h"

/**
* Spectrum smoother struct.
*/
typedef struct
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

	float release_coeff; //reference smoothing value
} SpectralSmoother;

/**
* Spectral smoothing proposed in 'Spectral subtraction with adaptive averaging of
* the gain function' but is not used yet.
*/
void spectrum_adaptive_time_smoothing(SpectralSmoother *self, float *prev_beta, float coeff)
{
	int k;
	float discrepancy, numerator = 0.f, denominator = 0.f;
	float beta_ts;
	float beta_smooth;
	float gamma_ts;

	for (k = 0; k <= self->half_fft_size; k++)
	{
		//These has to be magnitude spectrums
		numerator += fabs(self->signal_spectrum[k] - self->noise_spectrum[k]);
		denominator += self->noise_spectrum[k];
	}
	//this is the discrepancy of the spectum
	discrepancy = numerator / denominator;
	//beta is the adaptive coefficient
	beta_ts = MIN(discrepancy, 1.f);

	//Gamma is the smoothing coefficient of the adaptive factor beta
	if (*prev_beta < beta_ts)
	{
		gamma_ts = 0.f;
	}
	else
	{
		gamma_ts = coeff;
	}

	//Smoothing beta
	beta_smooth = gamma_ts * *(prev_beta) + (1.f - gamma_ts) * beta_ts;

	//copy current value to previous
	*prev_beta = beta_smooth;

	//Apply the adaptive smoothed beta over the signal
	for (k = 0; k <= self->half_fft_size; k++)
	{
		self->smoothed_spectrum[k] = (1.f - beta_smooth) * self->smoothed_spectrum_prev[k] + beta_smooth * self->smoothed_spectrum[k];
	}
}

/*
* Exponential decay coefficients for envelopes and adaptive noise profiling
* These must take into account the hop size as explained in the following paper
* FFT-BASED DYNAMIC RANGE COMPRESSION
*/
void get_release_coeff(SpectralSmoother *self, float release)
{
	if (release != 0.f) //This allows to turn off smoothing with 0 ms in order to use masking only
	{
		self->release_coeff = expf(-1000.f / (((release)*self->samp_rate) / self->hop));
	}
	else
	{
		self->release_coeff = 0.f; //This avoids incorrect results when moving sliders rapidly
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
			self->smoothed_spectrum[k] = self->release_coeff * self->smoothed_spectrum_prev[k] + (1.f - self->release_coeff) * self->smoothed_spectrum[k];
		}
	}
}

void s_s_run(SpectralSmoother *self, float release)
{
	get_release_coeff(self, release);

	memcpy(self->smoothed_spectrum, self->signal_spectrum, sizeof(float) * (self->half_fft_size + 1));

	apply_time_envelope(self);

	memcpy(self->smoothed_spectrum_prev, self->smoothed_spectrum, sizeof(float) * (self->half_fft_size + 1));
}

/**
* Reset dynamic arrays to zero.
*/
void s_s_reset(SpectralSmoother *self)
{
	//Reset all arrays
	initialize_array(self->signal_spectrum, 0.f, self->half_fft_size + 1);
	initialize_array(self->noise_spectrum, 0.f, self->half_fft_size + 1);
	initialize_array(self->smoothed_spectrum, 0.f, self->half_fft_size + 1);
	initialize_array(self->smoothed_spectrum_prev, 0.f, self->half_fft_size + 1);

	self->release_coeff = 0.f;
}

/**
* Gain estimator initialization and configuration.
*/
SpectralSmoother *
s_s_init(int fft_size, int samp_rate, int hop)
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
	s_s_reset(self);

	return self;
}

/**
* Free allocated memory.
*/
void s_s_free(SpectralSmoother *self)
{
	free(self->noise_spectrum);
	free(self->signal_spectrum);
	free(self->smoothed_spectrum);
	free(self->smoothed_spectrum_prev);
	free(self);
}