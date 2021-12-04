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

#include "fft_denoiser.h"
#include "gain_estimator.h"
#include "noise_estimator.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define WHITENING_DECAY_RATE 1000.f
#define WHITENING_FLOOR 0.02f

struct FFTDenoiser
{
	int fft_size;
	int half_fft_size;
	int samp_rate;
	int hop;

	float *fft_spectrum;
	float *processed_fft_spectrum;

	float tau;
	float wet_dry_target;
	float wet_dry;

	float *gain_spectrum;
	float *residual_spectrum;
	float *denoised_spectrum;
	float *whitened_residual_spectrum;

	GainEstimator *gain_estimation;
	NoiseEstimator *noise_estimation;
	NoiseProfile *noise_profile;

	float *residual_max_spectrum;
	float max_decay_rate;
	int whitening_window_count;
};

bool is_empty(float *spectrum, int N)
{
	int k;
	for (k = 0; k <= N; k++)
	{
		if (spectrum[k] > FLT_MIN)
		{
			return false;
		}
	}
	return true;
}

void fft_denoiser_update_wetdry_target(FFTDenoiser *self, bool enable)
{

	if (enable)
	{
		self->wet_dry_target = 1.f;
	}
	else
	{
		self->wet_dry_target = 0.f;
	}

	self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry) + FLT_MIN;
}

void fft_denoiser_soft_bypass(FFTDenoiser *self)
{
	int k;

	for (k = 0; k < self->fft_size; k++)
	{
		self->processed_fft_spectrum[k] = (1.f - self->wet_dry) * self->fft_spectrum[k] + self->processed_fft_spectrum[k] * self->wet_dry;
	}
}

void residual_spectrum_whitening(FFTDenoiser *self, float whitening_factor)
{
	self->whitening_window_count++;

	for (int k = 0; k < self->fft_size; k++)
	{
		if (self->whitening_window_count > 1.f)
		{
			self->residual_max_spectrum[k] = fmaxf(fmaxf(self->residual_spectrum[k], WHITENING_FLOOR), self->residual_max_spectrum[k] * self->max_decay_rate);
		}
		else
		{
			self->residual_max_spectrum[k] = fmaxf(self->residual_spectrum[k], WHITENING_FLOOR);
		}
	}

	for (int k = 0; k < self->fft_size; k++)
	{
		if (self->residual_spectrum[k] > FLT_MIN)
		{
			self->whitened_residual_spectrum[k] = self->residual_spectrum[k] / self->residual_max_spectrum[k];

			self->residual_spectrum[k] = (1.f - whitening_factor) * self->residual_spectrum[k] + whitening_factor * self->whitened_residual_spectrum[k];
		}
	}
}

void get_denoised_spectrum(FFTDenoiser *self)
{
	int k;

	for (k = 0; k < self->fft_size; k++)
	{
		self->denoised_spectrum[k] = self->fft_spectrum[k] * self->gain_spectrum[k];
	}
}

void get_residual_spectrum(FFTDenoiser *self, float whitening_factor)
{
	int k;

	for (k = 0; k < self->fft_size; k++)
	{
		self->residual_spectrum[k] = self->fft_spectrum[k] - self->denoised_spectrum[k];
	}

	if (whitening_factor > 0.f)
	{
		residual_spectrum_whitening(self, whitening_factor);
	}
}

void get_final_spectrum(FFTDenoiser *self, bool residual_listen, float reduction_amount)
{
	int k;

	if (residual_listen)
	{
		for (k = 0; k < self->fft_size; k++)
		{
			self->processed_fft_spectrum[k] = self->residual_spectrum[k];
		}
	}
	else
	{
		for (k = 0; k < self->fft_size; k++)
		{
			self->processed_fft_spectrum[k] = self->denoised_spectrum[k] +
											  self->residual_spectrum[k] * reduction_amount;
		}
	}
}

void fft_denoiser_run(FFTDenoiser *self, NoiseProfile *noise_profile, float *fft_spectrum, int enable, bool learn_noise, float whitening_factor,
					  float reduction_amount, bool residual_listen, float transient_threshold,
					  float masking_ceiling_limit, float release, float noise_rescale)
{
	fft_denoiser_update_wetdry_target(self, enable);

	memcpy(self->fft_spectrum, fft_spectrum, sizeof(float) * self->fft_size);

	if (!is_empty(self->fft_spectrum, self->half_fft_size))
	{
		if (learn_noise)
		{
			noise_estimation_run(self->noise_estimation, noise_profile, self->fft_spectrum);
		}
		else
		{
			if (is_noise_estimation_available(self->noise_estimation))
			{
				gain_estimation_run(self->gain_estimation, self->fft_spectrum, get_noise_profile(self->noise_profile), self->gain_spectrum, transient_threshold,
									masking_ceiling_limit, release, noise_rescale);

				get_denoised_spectrum(self);

				get_residual_spectrum(self, whitening_factor);

				get_final_spectrum(self, residual_listen, reduction_amount);
			}
		}
	}

	fft_denoiser_soft_bypass(self);

	memcpy(fft_spectrum, self->processed_fft_spectrum, sizeof(float) * self->fft_size);
}

void fft_denoiser_reset(FFTDenoiser *self)
{
	memset(self->fft_spectrum, 0.f, self->fft_size);
	memset(self->processed_fft_spectrum, 0.f, self->fft_size);
	memset(self->gain_spectrum, 1.f, self->fft_size);

	memset(self->residual_max_spectrum, 0.f, self->fft_size);
	memset(self->denoised_spectrum, 0.f, self->fft_size);
	memset(self->residual_spectrum, 0.f, self->fft_size);
	memset(self->whitened_residual_spectrum, 0.f, self->fft_size);
	memset(self->gain_spectrum, 0.f, self->fft_size);

	self->whitening_window_count = 0.f;
}

FFTDenoiser *fft_denoiser_initialize(int samp_rate, int fft_size, int overlap_factor)
{
	FFTDenoiser *self = (FFTDenoiser *)malloc(sizeof(FFTDenoiser));

	self->fft_size = fft_size;
	self->half_fft_size = self->fft_size / 2;
	self->hop = self->fft_size / overlap_factor;
	self->samp_rate = samp_rate;

	self->fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->processed_fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));

	self->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->samp_rate));
	self->wet_dry = 0.f;

	self->residual_max_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->max_decay_rate = expf(-1000.f / (((WHITENING_DECAY_RATE)*self->samp_rate) / self->hop));

	self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->gain_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->whitened_residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));

	fft_denoiser_reset(self);

	self->noise_estimation = noise_estimation_initialize(self->fft_size);

	return self;
}

void fft_denoiser_free(FFTDenoiser *self)
{
	free(self->fft_spectrum);
	free(self->processed_fft_spectrum);
	free(self->gain_spectrum);
	free(self->residual_spectrum);
	free(self->whitened_residual_spectrum);
	free(self->denoised_spectrum);
	free(self->residual_max_spectrum);
	noise_estimation_free(self->noise_estimation);
	free(self);
}
