/*
libspecbleach - A spectral processing library

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

#include "stft_buffer.h"
#include <stdlib.h>
#include <string.h>

struct StftBuffer {
  uint32_t read_position;
  uint32_t start_position;
  uint32_t stft_frame_size;
  uint32_t block_step;

  // TODO (luciano/todo): replace FIFO buffers with one single lock free Queue
  float *in_fifo;
  float *out_fifo;
};

StftBuffer *stft_buffer_initialize(const uint32_t stft_frame_size,
                                   const uint32_t start_position,
                                   const uint32_t block_step) {
  StftBuffer *self = (StftBuffer *)calloc(1U, sizeof(StftBuffer));

  self->stft_frame_size = stft_frame_size;
  self->start_position = start_position;
  self->block_step = block_step;
  self->read_position = self->start_position;
  self->in_fifo = (float *)calloc(self->stft_frame_size, sizeof(float));
  self->out_fifo = (float *)calloc(self->stft_frame_size, sizeof(float));

  return self;
}

void stft_buffer_free(StftBuffer *self) {
  free(self->in_fifo);
  free(self->out_fifo);
  free(self);
}

bool stft_buffer_fill(StftBuffer *self, const float input_sample,
                      float *output_sample) {
  if (!output_sample) {
    return false;
  }

  self->in_fifo[self->read_position] = input_sample;
  *output_sample = self->out_fifo[self->read_position - self->start_position];
  self->read_position++;

  // Is it full?
  if (self->read_position == self->stft_frame_size) {
    return true;
  }

  return false;
}

bool stft_buffer_advance_block(StftBuffer *self,
                               const float *reconstructed_signal) {
  if (!reconstructed_signal) {
    return false;
  }

  self->read_position = self->start_position;

  memcpy(self->in_fifo, &self->in_fifo[self->block_step],
         sizeof(float) * self->start_position);

  memcpy(self->out_fifo, reconstructed_signal,
         sizeof(float) * self->block_step);

  return true;
}

float *get_full_buffer_block(StftBuffer *self) { return self->in_fifo; }