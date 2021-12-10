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

STFTProcessor *stft_processor_initialize(FFTDenoiser *fft_denoiser, int fft_size, int overlap_factor);
void stft_processor_free(STFTProcessor *self);
void set_spectral_size(STFTProcessor *self, int fft_size);
int stft_processor_get_latency(STFTProcessor *self);
void stft_processor_run(STFTProcessor *self, NoiseProfile *noise_profile, int n_samples, const float *input, float *output);

#endif