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

typedef enum ZeroPaddingType {
  NEXT_POWER_OF_TWO = 0,
  FIXED_AMOUNT = 1,
  NO_PADDING = 2,
} ZeroPaddingType;

typedef struct FftTransform FftTransform;

FftTransform *fft_transform_initialize(uint32_t sample_rate,
                                       float frame_size_ms,
                                       ZeroPaddingType padding_type);
FftTransform *fft_transform_initialize_bins(uint32_t fft_size);
void fft_transform_free(FftTransform *self);
bool fft_load_input_samples(FftTransform *self, const float *input);
bool fft_get_output_samples(FftTransform *self, float *output);
uint32_t get_fft_size(FftTransform *self);
uint32_t get_frame_size(FftTransform *self);
uint32_t get_fft_real_spectrum_size(FftTransform *self);
bool compute_forward_fft(FftTransform *self);
bool compute_backward_fft(FftTransform *self);
float *get_fft_input_buffer(FftTransform *self);
float *get_fft_output_buffer(FftTransform *self);

#endif