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
* \file fft_processor.c
* \author Luciano Dato
* \brief Contains an abstraction for a single fft spectrum denoising
*/

#define WHITENING_DECAY_RATE 1000.f //Deacay in ms for max spectrum for whitening
#define WHITENING_FLOOR 0.02f		//Minumum max value posible

#include "extra_functions.h"
#include "gain_estimator.c"
#include "noise_estimator.c"

/**
* FFT processor struct.
*/
typedef struct
{
	//General parameters
	int fft_size;
	int fft_size_2;
	int samp_rate;
	int hop;

	//Ensemble related
	float *fft_spectrum;
	float *processed_fft_spectrum;

	//Soft bypass
	float tau;			  //time constant for soft bypass
	float wet_dry_target; //softbypass target for softbypass
	float wet_dry;		  //softbypass coeff

	//Spectrum information arrays
	float *power_spectrum;
	float *phase_spectrum;
	float *magnitude_spectrum;

	float *gain_spectrum; //definitive reduction gain
	float *residual_spectrum;
	float *denoised_spectrum;
	float *whitened_residual_spectrum;

	GainEstimator *gain_estimation;
	NoiseEstimator *noise_estimation;

	//whitening related
	float *residual_max_spectrum;
	float max_decay_rate;
	float whitening_window_count;
} FFTProcessor;

/**
* Updates the wet/dry mixing coefficient.
*/
void fft_processor_update_wetdry_target(FFTProcessor *self, bool enable)
{
	//Softbypass targets in case of disabled or enabled
	if (enable)
	{ //if enabled
		self->wet_dry_target = 1.f;
	}
	else
	{ //if disabled
		self->wet_dry_target = 0.f;
	}
	//Interpolate parameters over time softly to bypass without clicks or pops
	self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry) + FLT_MIN;
}

/**
* Mixes unprocessed and processed signal to bypass softly.
*/
void fft_processor_soft_bypass(FFTProcessor *self)
{
	int k;

	for (k = 0; k < self->fft_size; k++)
	{
		self->processed_fft_spectrum[k] = (1.f - self->wet_dry) * self->fft_spectrum[k] + self->processed_fft_spectrum[k] * self->wet_dry;
	}
}

/**
* Whitens the spectrum adaptively as proposed in 'Adaptive whitening for improved
* real-time audio onset detection' by Stowell and Plumbley. The idea here is that when
* residual noise resembles white noise the ear is able to precieve it as not so annoying.
* It uses a temporal max value for each bin and a decay factor as the memory regulator of
* that maximun value.
*/
void residual_spectrum_whitening(FFTProcessor *self, float whitening_factor)
{
	self->whitening_window_count++;

	for (int k = 0; k < self->fft_size; k++)
	{
		if (self->whitening_window_count > 1.f)
		{
			self->residual_max_spectrum[k] = MAX(MAX(self->residual_spectrum[k], WHITENING_FLOOR), self->residual_max_spectrum[k] * self->max_decay_rate);
		}
		else
		{
			self->residual_max_spectrum[k] = MAX(self->residual_spectrum[k], WHITENING_FLOOR);
		}
	}

	for (int k = 0; k < self->fft_size; k++)
	{
		if (self->residual_spectrum[k] > FLT_MIN)
		{
			//Get whitened spectrum
			self->whitened_residual_spectrum[k] = self->residual_spectrum[k] / self->residual_max_spectrum[k];

			//Interpolate between whitened and non whitened residual
			self->residual_spectrum[k] = (1.f - whitening_factor) * self->residual_spectrum[k] + whitening_factor * self->whitened_residual_spectrum[k];
		}
	}
}

/**
* Applies the filter to the complex spectrum and gets the clean signal.
*/
void get_denoised_spectrum(FFTProcessor *self)
{
	int k;

	//Apply the computed gain to the signal and store it in denoised array
	for (k = 0; k < self->fft_size; k++)
	{
		self->denoised_spectrum[k] = self->fft_spectrum[k] * self->gain_spectrum[k];
	}
}

/**
* Gets the residual signal of the reduction.
*/
void get_residual_spectrum(FFTProcessor *self, float whitening_factor)
{
	int k;

	//Residual signal
	for (k = 0; k < self->fft_size; k++)
	{
		self->residual_spectrum[k] = self->fft_spectrum[k] - self->denoised_spectrum[k];
	}

	//Whitening (residual spectrum more similar to white noise)
	if (whitening_factor > 0.f)
	{
		residual_spectrum_whitening(self, whitening_factor);
	}
}

/**
* Mixes the cleaned signal with the residual taking into account the reduction configured
* by the user. Outputs the final signal or the residual only.
*/
void get_final_spectrum(FFTProcessor *self, bool residual_listen, float reduction_amount)
{
	int k;

	//OUTPUT RESULTS using smooth bypass and parametric subtraction
	if (residual_listen)
	{
		//Output noise only
		for (k = 0; k < self->fft_size; k++)
		{
			self->processed_fft_spectrum[k] = self->residual_spectrum[k];
		}
	}
	else
	{
		//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k < self->fft_size; k++)
		{
			self->processed_fft_spectrum[k] = self->denoised_spectrum[k] +
											  self->residual_spectrum[k] * reduction_amount;
		}
	}
}

/**
* Runs the fft processing for current block.
*/
void fft_processor_run(FFTProcessor *self, float *fft_spectrum, int enable, bool learn_noise, float whitening_factor,
					   float reduction_amount, bool residual_listen, float transient_threshold,
					   float masking_ceiling_limit, float release, float noise_rescale)
{
	fft_processor_update_wetdry_target(self, enable);

	memcpy(self->fft_spectrum, fft_spectrum, sizeof(float) * self->fft_size);

	//First get the power, magnitude and phase spectrum
	get_info_from_bins(self->power_spectrum, self->magnitude_spectrum,
					   self->phase_spectrum, self->fft_size_2,
					   self->fft_size, self->fft_spectrum);

	//If the spectrum is not silence
	if (!is_empty(self->power_spectrum, self->fft_size_2))
	{
		/*If selected estimate noise spectrum is based on selected portion of signal
        *do not process the signal
        */
		if (learn_noise)
		{
			//LEARN NOISE (Using power spectrum)
			noise_estimation_run(self->noise_estimation, self->power_spectrum);
		}
		else
		{
			//REDUCE NOISE OR LISTEN TO THE RESIDUAL
			if (is_noise_estimation_available(self->noise_estimation))
			{
				gain_estimation_run(self->gain_estimation, self->power_spectrum, self->gain_spectrum, transient_threshold,
									masking_ceiling_limit, release, noise_rescale);

				get_denoised_spectrum(self);

				get_residual_spectrum(self, whitening_factor);

				get_final_spectrum(self, residual_listen, reduction_amount);
			}
		}
	}

	//If bypassed mix unprocessed and processed signal softly
	fft_processor_soft_bypass(self);

	//Copy the processed spectrum to fft_spectrum
	memcpy(fft_spectrum, self->processed_fft_spectrum, sizeof(float) * self->fft_size);
}

/**
* Reset dynamic arrays to zero.
*/
void fft_processor_reset(FFTProcessor *self)
{
	//Reset all arrays
	initialize_array(self->fft_spectrum, 0.f, self->fft_size);
	initialize_array(self->processed_fft_spectrum, 0.f, self->fft_size);
	initialize_array(self->gain_spectrum, 1.f, self->fft_size);

	initialize_array(self->power_spectrum, 0.f, self->fft_size_2 + 1);
	initialize_array(self->magnitude_spectrum, 0.f, self->fft_size_2 + 1);
	initialize_array(self->phase_spectrum, 0.f, self->fft_size_2 + 1);

	initialize_array(self->residual_max_spectrum, 0.f, self->fft_size);
	initialize_array(self->denoised_spectrum, 0.f, self->fft_size);
	initialize_array(self->residual_spectrum, 0.f, self->fft_size);
	initialize_array(self->whitened_residual_spectrum, 0.f, self->fft_size);
	initialize_array(self->gain_spectrum, 0.f, self->fft_size);

	self->whitening_window_count = 0.f;
}

/**
* FFT processor initialization and configuration.
*/
FFTProcessor *
fft_processor_initialize(int fft_size, int samp_rate, int hop)
{
	//Allocate object
	FFTProcessor *self = (FFTProcessor *)malloc(sizeof(FFTProcessor));

	//Configuration
	self->fft_size = fft_size;
	self->fft_size_2 = self->fft_size / 2;
	self->samp_rate = samp_rate;
	self->hop = hop;

	//spectrum allocation
	self->fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->processed_fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));

	//soft bypass
	self->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->samp_rate));
	self->wet_dry = 0.f;

	//whitening related
	self->residual_max_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->max_decay_rate = expf(-1000.f / (((WHITENING_DECAY_RATE)*self->samp_rate) / self->hop));

	//final ensemble related
	self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->gain_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->whitened_residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));

	//Arrays for getting bins info
	self->power_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));
	self->magnitude_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));
	self->phase_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));

	//Reset all values
	fft_processor_reset(self);

	//Noise estimator related
	self->noise_estimation = noise_estimation_initialize(self->fft_size);

	return self;
}

/**
* Free allocated memory.
*/
void fft_processor_free(FFTProcessor *self)
{
	free(self->fft_spectrum);
	free(self->processed_fft_spectrum);
	free(self->power_spectrum);
	free(self->magnitude_spectrum);
	free(self->phase_spectrum);
	free(self->gain_spectrum);
	free(self->residual_spectrum);
	free(self->whitened_residual_spectrum);
	free(self->denoised_spectrum);
	free(self->residual_max_spectrum);
	noise_estimation_free(self->noise_estimation);
	free(self);
}
