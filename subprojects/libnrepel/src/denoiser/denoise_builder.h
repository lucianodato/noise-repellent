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

#ifndef DENOISE_BUILDER_H
#define DENOISE_BUILDER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct DenoiseBuilder DenoiseBuilder;

DenoiseBuilder *denoise_builder_initialize(uint32_t half_fft_size,
                                           uint32_t fft_size, uint32_t hop,
                                           uint32_t sample_rate);
void denoise_builder_free(DenoiseBuilder *self);
void denoise_build(DenoiseBuilder *self, bool residual_listening,
                   float reduction_amout, float whitening_factor,
                   const float *gain_spectrum, float *fft_spectrum);

#endif