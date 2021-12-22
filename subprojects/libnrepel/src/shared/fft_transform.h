/*
noise-repellent -- Noise Reduction LV2

Copyright 2021 Luciano Dato <lucianodato@gmail.com>

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

#ifndef FFT_TRANSFORM_H
#define FFT_TRANSFORM_H

#include <stdbool.h>
#include <stdint.h>

typedef struct FftTransform FftTransform;

FftTransform *fft_transform_initialize();
void fft_transform_free(FftTransform *self);
bool load_input_samples(FftTransform *self, const float *input);
uint32_t get_fft_size(FftTransform *self);
uint32_t get_real_spectrum_size(FftTransform *self);
bool compute_forward_fft(FftTransform *self);
bool compute_backward_fft(FftTransform *self);
float *get_fft_input_buffer(FftTransform *self);
float *get_fft_output_buffer(FftTransform *self);

#endif