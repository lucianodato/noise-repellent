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
* \file transient_detector.c
* \author Luciano Dato
* \brief Contains a transient detector abstraction
*/

#include "extra_functions.h"
#include <stdbool.h>

#define TP_UPPER_LIMIT 5.f //This correspond to the upper limit of the adaptive threshold multiplier. Should be the same as the ttl configured one

/**
* Gain estimation struct.
*/
typedef struct
{
	//General parameters
	int fft_size;
	int half_fft_size;

	float *spectrum;

	//Transient preservation related
	float *transient_preserv_prev; //previous frame smoothed power spectrum for envelopes
	float tp_r_mean;
	bool transient_present;
	float tp_window_count;
} TransientDetector;

/**
* Outputs the spectral flux between two spectrums.
* \param spectrum the current power spectrum
* \param spectrum_prev the previous power spectrum
* \param N the size of the spectrum (half the fft size plus 1)
*/
float spectral_flux(float *spectrum, float *spectrum_prev, float N)
{
	int i;
	float spectral_flux = 0.f;
	float temp;

	for (i = 0; i <= N; i++)
	{
		temp = sqrtf(spectrum[i]) - sqrtf(spectrum_prev[i]); //Recieves power spectrum uses magnitude
		spectral_flux += (temp + fabs(temp)) / 2.f;
	}
	return spectral_flux;
}

/**
* Outputs the high frequency content of the spectrum.
* \param spectrum the current power spectrum
* \param N the size of the spectrum (half the fft size plus 1)
*/
float high_frequency_content(float *spectrum, float N)
{
	int i;
	float sum = 0.f;

	for (i = 0; i <= N; i++)
	{
		sum += i * spectrum[i];
	}
	return sum / (float)(N + 1);
}

/**
* Transient detection using a rolling mean thresholding over the spectral flux of
* the signal. Using more heuristics like high frequency content and others like the ones
* anylised by Dixon in 'Simple Spectrum-Based Onset Detection' would be better. Onset
* detection is explained thoroughly in 'A tutorial on onset detection in music signals' * by Bello.
*/
bool t_d_run(TransientDetector *self, float transient_threshold)
{
	float adapted_threshold, reduction_function;

	//Transient protection by forcing wiener filtering when an onset is detected
	reduction_function = spectral_flux(self->spectrum, self->transient_preserv_prev, self->half_fft_size);
	//reduction_function = high_frequency_content(self->spectrum, self->half_fft_size);

	self->tp_window_count += 1.f;

	if (self->tp_window_count > 1.f)
	{
		self->tp_r_mean += ((reduction_function - self->tp_r_mean) / self->tp_window_count);
	}
	else
	{
		self->tp_r_mean = reduction_function;
	}

	adapted_threshold = (TP_UPPER_LIMIT - transient_threshold) * self->tp_r_mean;

	memcpy(self->transient_preserv_prev, self->spectrum, sizeof(float) * (self->half_fft_size + 1));

	if (reduction_function > adapted_threshold)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/**
* Reset dynamic arrays to zero.
*/
void t_d_reset(TransientDetector *self)
{
	//Reset all arrays
	initialize_array(self->spectrum, 0.f, self->half_fft_size + 1);
	initialize_array(self->transient_preserv_prev, 0.f, self->half_fft_size + 1);

	self->tp_window_count = 0.f;
	self->tp_r_mean = 0.f;
	self->transient_present = false;
}

/**
* Masking estimator initialization and configuration.
*/
TransientDetector *
t_d_init(int fft_size)
{
	//Allocate object
	TransientDetector *self = (TransientDetector *)malloc(sizeof(TransientDetector));

	//Configuration
	self->fft_size = fft_size;
	self->half_fft_size = self->fft_size / 2;

	//spectrum allocation
	self->spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->transient_preserv_prev = (float *)calloc((self->half_fft_size + 1), sizeof(float));

	//Reset all values
	t_d_reset(self);

	return self;
}

/**
* Free allocated memory.
*/
void t_d_free(TransientDetector *self)
{
	free(self->spectrum);
	free(self->transient_preserv_prev);
	free(self);
}