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

#include "fft_transform.h"
#include "../shared/configurations.h"
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>

struct FftTransform {
  fftwf_plan forward;
  fftwf_plan backward;

  uint32_t fft_size;
  float *input_fft_buffer;
  float *output_fft_buffer;
};

FftTransform *fft_transform_initialize() {
  FftTransform *self = (FftTransform *)calloc(1U, sizeof(FftTransform));

  self->fft_size = FFT_SIZE;

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

void fft_transform_free(FftTransform *self) {
  free(self->input_fft_buffer);
  free(self->output_fft_buffer);
  fftwf_destroy_plan(self->forward);
  fftwf_destroy_plan(self->backward);

  free(self);
}

uint32_t get_fft_size(FftTransform *self) { return self->fft_size; }
uint32_t get_real_spectrum_size(FftTransform *self) {
  return self->fft_size / 2U + 1U;
}

bool fft_load_input_samples(FftTransform *self, const float *input) {
  if (!self || !input) {
    return false;
  }

  memcpy(self->input_fft_buffer, input, sizeof(float) * self->fft_size);

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