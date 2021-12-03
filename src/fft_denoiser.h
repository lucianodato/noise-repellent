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
* \file fft_denoiser.h
* \author Luciano Dato
* \brief Contains an abstraction for a single fft spectrum denoising
*/

#ifndef FFT_DENOISER_H
#define FFT_DENOISER_H

#include <stdbool.h>

typedef struct FFTDenoiser FFTDenoiser;

bool is_empty(float *spectrum, int N);
void fft_denoiser_update_wetdry_target(FFTDenoiser *self, bool enable);
void fft_denoiser_soft_bypass(FFTDenoiser *self);
void residual_spectrum_whitening(FFTDenoiser *self, float whitening_factor);
void get_denoised_spectrum(FFTDenoiser *self);
void get_residual_spectrum(FFTDenoiser *self, float whitening_factor);
void get_final_spectrum(FFTDenoiser *self, bool residual_listen, float reduction_amount);
void fft_denoiser_run(FFTDenoiser *self, float *fft_spectrum, int enable, bool learn_noise, float whitening_factor,
					  float reduction_amount, bool residual_listen, float transient_threshold,
					  float masking_ceiling_limit, float release, float noise_rescale);
void fft_denoiser_reset(FFTDenoiser *self);
FFTDenoiser *fft_denoiser_initialize(int fft_size, int half_fft_size, int samp_rate, int hop);
void fft_denoiser_free(FFTDenoiser *self);

#endif