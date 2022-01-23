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

#include "fft_transform.h"
#include "../shared/configurations.h"
#include "../shared/general_utils.h"

#include <fftw3.h>
#include <stdlib.h>
#include <string.h>

static uint32_t calculate_fft_size(FftTransform *self);

struct FftTransform {
  fftwf_plan forward;
  fftwf_plan backward;

  uint32_t fft_size;
  uint32_t frame_size;
  ZeroPaddingType padding_type;
  uint32_t padding_amount;
  float *input_fft_buffer;
  float *output_fft_buffer;
};

FftTransform *fft_transform_initialize(const uint32_t sample_rate,
                                       const float frame_size_ms,
                                       const ZeroPaddingType padding_type) {
  FftTransform *self = (FftTransform *)calloc(1U, sizeof(FftTransform));

  self->padding_type = padding_type;
  self->frame_size = (uint32_t)((frame_size_ms / 1000.F) * (float)sample_rate);

  self->fft_size = calculate_fft_size(self);

  self->input_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->output_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->forward =
      fftwf_plan_r2r_1d((int)self->fft_size, self->input_fft_buffer,
                        self->output_fft_buffer, FFTW_FORWARD, FFTW_ESTIMATE);
  self->backward =
      fftwf_plan_r2r_1d((int)self->fft_size, self->output_fft_buffer,
                        self->input_fft_buffer, FFTW_BACKWARD, FFTW_ESTIMATE);

  return self;
}

FftTransform *fft_transform_initialize_bins(const uint32_t fft_size) {
  FftTransform *self = (FftTransform *)calloc(1U, sizeof(FftTransform));

  self->fft_size = fft_size;
  self->frame_size = self->fft_size;

  self->input_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->output_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->forward =
      fftwf_plan_r2r_1d((int)self->fft_size, self->input_fft_buffer,
                        self->output_fft_buffer, FFTW_FORWARD, FFTW_ESTIMATE);
  self->backward =
      fftwf_plan_r2r_1d((int)self->fft_size, self->output_fft_buffer,
                        self->input_fft_buffer, FFTW_BACKWARD, FFTW_ESTIMATE);

  return self;
}

static uint32_t calculate_fft_size(FftTransform *self) {
  switch (self->padding_type) {
  case NO_PADDING: {
    self->padding_amount = 0;
    return get_next_divisible_two((int)self->frame_size);
  }
  case NEXT_POWER_OF_TWO: {
    uint32_t next_power_of_two = get_next_power_two((int)self->frame_size);
    self->padding_amount = next_power_of_two - self->frame_size;
    return next_power_of_two;
  }
  case FIXED_AMOUNT: {
    self->padding_amount = ZEROPADDING_AMOUNT;
    return get_next_divisible_two(
        (int)(self->frame_size + self->padding_amount));
  }
  default:
    return get_next_divisible_two((int)self->frame_size);
  }
}

void fft_transform_free(FftTransform *self) {
  free(self->input_fft_buffer);
  free(self->output_fft_buffer);
  fftwf_destroy_plan(self->forward);
  fftwf_destroy_plan(self->backward);

  free(self);
}

uint32_t get_fft_size(FftTransform *self) { return self->fft_size; }
uint32_t get_frame_size(FftTransform *self) { return self->frame_size; }
uint32_t get_fft_real_spectrum_size(FftTransform *self) {
  return self->fft_size / 2U + 1U;
}

bool fft_load_input_samples(FftTransform *self, const float *input) {
  if (!self || !input) {
    return false;
  }

  memcpy(&self->input_fft_buffer[self->padding_amount / 2U], input,
         sizeof(float) * self->frame_size);

  return true;
}

bool fft_get_output_samples(FftTransform *self, float *output) {
  if (!self || !output) {
    return false;
  }

  memcpy(output, &self->input_fft_buffer[self->padding_amount / 2U],
         sizeof(float) * self->frame_size);

  return true;
}

bool compute_forward_fft(FftTransform *self) {
  if (!self) {
    return false;
  }

  fftwf_execute(self->forward);

  return true;
}

bool compute_backward_fft(FftTransform *self) {
  if (!self) {
    return false;
  }

  fftwf_execute(self->backward);

  return true;
}

float *get_fft_input_buffer(FftTransform *self) {
  return self->input_fft_buffer;
}

float *get_fft_output_buffer(FftTransform *self) {
  return self->output_fft_buffer;
}