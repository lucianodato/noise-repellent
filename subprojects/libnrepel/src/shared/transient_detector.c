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

#include "transient_detector.h"
#include "configurations.h"
#include "spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct TransientDetector {
  uint32_t fft_size;
  uint32_t real_spectrum_size;

  float *previous_spectrum;
  float rolling_mean;
  bool transient_present;
  uint32_t window_count;
};

TransientDetector *transient_detector_initialize(const uint32_t fft_size) {
  TransientDetector *self =
      (TransientDetector *)calloc(1U, sizeof(TransientDetector));

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;

  self->previous_spectrum =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

  self->window_count = 0U;
  self->rolling_mean = 0.F;
  self->transient_present = false;

  return self;
}

void transient_detector_free(TransientDetector *self) {
  free(self->previous_spectrum);
  free(self);
}

bool transient_detector_run(TransientDetector *self,
                            const float transient_threshold,
                            const float *spectrum) {
  const float reduction_function = spectral_flux(
      spectrum, self->previous_spectrum, self->real_spectrum_size);

  self->window_count += 1U;

  if (self->window_count > 1U) {
    self->rolling_mean +=
        ((reduction_function - self->rolling_mean) / (float)self->window_count);
  } else {
    self->rolling_mean = reduction_function;
  }

  const float adapted_threshold =
      (UPPER_LIMIT - transient_threshold) * self->rolling_mean;

  memcpy(self->previous_spectrum, spectrum,
         sizeof(float) * self->real_spectrum_size);

  if (reduction_function > adapted_threshold) {
    return true;
  }
  return false;
}