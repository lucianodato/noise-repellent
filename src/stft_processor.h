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

#ifndef STFT_PROCESSOR_H
#define STFT_PROCESSOR_H

#include "fft_denoiser.h"
#include "noise_profile.h"
#include <stdbool.h>

typedef struct STFTProcessor STFTProcessor;

void set_spectral_size(STFTProcessor *self, int fft_size);
void fft_window(float *window, int N, int window_type);
void stft_processor_pre_and_post_window(STFTProcessor *self);
void stft_processor_analysis(STFTProcessor *self);
void stft_processor_synthesis(STFTProcessor *self);
int stft_processor_get_latency(STFTProcessor *self);
void stft_processor_run(STFTProcessor *self, NoiseProfile *noise_profile, int n_samples, const float *input, float *output,
						int enable, int learn_noise, float whitening_factor, float reduction_amount,
						bool residual_listen, float transient_threshold, float masking_ceiling_limit,
						float release, float noise_rescale);
void stft_processor_reset(STFTProcessor *self);
STFTProcessor *stft_processor_initialize(FFTDenoiser *fft_denoiser, int sample_rate, int fft_size, int overlap_factor);
void stft_processor_free(STFTProcessor *self);
void get_info_from_bins(float *fft_p2, float *fft_magnitude, float *fft_phase,
						int fft_size_2, int fft_size, float *fft_buffer);

#endif