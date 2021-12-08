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

#ifndef FFT_DENOISER_H
#define FFT_DENOISER_H

#include "denoise_parameters.h"
#include "noise_profile.h"
#include <float.h>
#include <stdbool.h>

typedef struct FFTDenoiser FFTDenoiser;

bool is_empty(float *spectrum, int N);
void fft_denoiser_update_wetdry_target(FFTDenoiser *self, bool enable);
void fft_denoiser_soft_bypass(FFTDenoiser *self);
void residual_spectrum_whitening(FFTDenoiser *self, float whitening_factor);
void get_denoised_spectrum(FFTDenoiser *self);
void get_residual_spectrum(FFTDenoiser *self, float whitening_factor);
void get_final_spectrum(FFTDenoiser *self, bool residual_listen, float reduction_amount);
void fft_denoiser_run(FFTDenoiser *self, NoiseProfile *noise_profile, float *fft_spectrum);
FFTDenoiser *fft_denoiser_initialize(int sample_rate, int fft_size, int overlap_factor);
void fft_denoiser_free(FFTDenoiser *self);
void load_denoise_parameters(FFTDenoiser *self, DenoiseParameters *new_parameters);

#endif