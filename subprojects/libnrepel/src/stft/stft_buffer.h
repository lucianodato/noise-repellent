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

#ifndef STFT_BUFFER_H
#define STFT_BUFFER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct StftBuffer StftBuffer;
StftBuffer *stft_buffer_initialize(uint32_t stft_frame_size,
                                   uint32_t start_position,
                                   uint32_t block_step);
void stft_buffer_free(StftBuffer *self);
bool stft_buffer_fill(StftBuffer *self, float input_sample,
                      float *output_sample);
bool stft_buffer_advance_block(StftBuffer *self,
                               const float *reconstructed_signal);
float *get_full_buffer_block(StftBuffer *self);

#endif