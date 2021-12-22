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

#include "stft_processor.h"
#include "../shared/configurations.h"
#include "../shared/fft_transform.h"
#include "../shared/spectral_features.h"
#include "stft_windows.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void stft_analysis(StftProcessor *self);
static void stft_synthesis(StftProcessor *self);
static void stft_overlap_add(StftProcessor *self);
static void stft_update_buffers(StftProcessor *self);

typedef struct {
  float overlap_scale_factor;

  uint32_t input_latency;
  uint32_t read_position;
  uint32_t hop;
  uint32_t overlap_factor;
  uint32_t buffer_size;

  float *in_fifo;
  float *out_fifo;
  float *output_accumulator;

} SamplesBuffer;

struct StftProcessor {
  FftTransform *fft_transform;
  SamplesBuffer stft_buffer;
  StftWindows *stft_windows;
};

StftProcessor *stft_processor_initialize() {
  StftProcessor *self = (StftProcessor *)calloc(1U, sizeof(StftProcessor));

  self->fft_transform = fft_transform_initialize();

  self->stft_buffer.buffer_size = get_fft_size(self->fft_transform);

  self->stft_buffer.overlap_factor = OVERLAP_FACTOR;
  self->stft_buffer.hop =
      self->stft_buffer.buffer_size / self->stft_buffer.overlap_factor;
  self->stft_buffer.input_latency =
      self->stft_buffer.buffer_size - self->stft_buffer.hop;
  self->stft_buffer.read_position = self->stft_buffer.input_latency;
  self->stft_buffer.in_fifo =
      (float *)calloc(self->stft_buffer.buffer_size, sizeof(float));
  self->stft_buffer.out_fifo =
      (float *)calloc(self->stft_buffer.buffer_size, sizeof(float));
  self->stft_buffer.output_accumulator =
      (float *)calloc(self->stft_buffer.buffer_size * 2U, sizeof(float));

  self->stft_windows = stft_window_initialize(self->stft_buffer.buffer_size,
                                              self->stft_buffer.overlap_factor);

  return self;
}

void stft_processor_free(StftProcessor *self) {
  free(self->stft_buffer.in_fifo);
  free(self->stft_buffer.out_fifo);
  free(self->stft_buffer.output_accumulator);

  fft_transform_free(self->fft_transform);

  stft_window_free(self->stft_windows);

  free(self);
}

uint32_t get_stft_latency(StftProcessor *self) {
  return self->stft_buffer.input_latency;
}
uint32_t get_buffer_size(StftProcessor *self) {
  return self->stft_buffer.buffer_size;
}
uint32_t get_overlap_factor(StftProcessor *self) {
  return self->stft_buffer.overlap_factor;
}
uint32_t get_spectral_processing_size(StftProcessor *self) {
  return get_real_spectrum_size(self->fft_transform);
}

bool stft_processor_run(StftProcessor *self,
                        spectral_processing *spectral_processing,
                        void *spectral_processor,
                        const uint32_t number_of_samples, const float *input,
                        float *output) {
  if (!self || !spectral_processing || !spectral_processor || !input ||
      !output || number_of_samples <= 0U) {
    return false;
  }

  for (uint32_t k = 0U; k < number_of_samples; k++) {
    self->stft_buffer.in_fifo[self->stft_buffer.read_position] = input[k];
    output[k] = self->stft_buffer.out_fifo[self->stft_buffer.read_position -
                                           self->stft_buffer.input_latency];
    self->stft_buffer.read_position++;

    if (self->stft_buffer.read_position >= self->stft_buffer.buffer_size) {
      self->stft_buffer.read_position = self->stft_buffer.input_latency;

      load_input_samples(self->fft_transform, self->stft_buffer.in_fifo);

      stft_analysis(self);

      spectral_processing(spectral_processor,
                          get_fft_output_buffer(self->fft_transform));

      stft_synthesis(self);

      stft_overlap_add(self);

      stft_update_buffers(self);
    }
  }

  return true;
}

static void stft_analysis(StftProcessor *self) {
  apply_window(self->stft_windows, get_fft_input_buffer(self->fft_transform),
               self->stft_buffer.buffer_size, INPUT_WINDOW);

  compute_forward_fft(self->fft_transform);
}

static void stft_synthesis(StftProcessor *self) {
  compute_backward_fft(self->fft_transform);

  apply_window(self->stft_windows, get_fft_input_buffer(self->fft_transform),
               self->stft_buffer.buffer_size, OUTPUT_WINDOW);
}

static void stft_overlap_add(StftProcessor *self) {

  for (uint32_t k = 0U; k < self->stft_buffer.buffer_size; k++) {
    self->stft_buffer.output_accumulator[k] +=
        get_fft_input_buffer(self->fft_transform)[k];
  }

  memcpy(self->stft_buffer.out_fifo, self->stft_buffer.output_accumulator,
         sizeof(float) * self->stft_buffer.hop);
}

static void stft_update_buffers(StftProcessor *self) {
  memmove(self->stft_buffer.output_accumulator,
          &self->stft_buffer.output_accumulator[self->stft_buffer.hop],
          self->stft_buffer.buffer_size * sizeof(float));

  memcpy(self->stft_buffer.in_fifo,
         &self->stft_buffer.in_fifo[self->stft_buffer.hop],
         sizeof(float) * self->stft_buffer.input_latency);
}