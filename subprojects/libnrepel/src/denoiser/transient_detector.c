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

#include "transient_detector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define TP_UPPER_LIMIT 5.f

static float spectral_flux(const float *spectrum, float *spectrum_prev,
                           uint32_t half_fft_size);

struct TransientDetector {
  uint32_t fft_size;
  uint32_t half_fft_size;

  float *previous_spectrum;
  float rolling_mean;
  bool transient_present;
  float window_count;
};

TransientDetector *transient_detector_initialize(const uint32_t fft_size) {
  TransientDetector *self =
      (TransientDetector *)calloc(1, sizeof(TransientDetector));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;

  self->previous_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->window_count = 0.f;
  self->rolling_mean = 0.f;
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
  const float reduction_function =
      spectral_flux(spectrum, self->previous_spectrum, self->half_fft_size);

  self->window_count += 1.f;

  if (self->window_count > 1.f) {
    self->rolling_mean +=
        ((reduction_function - self->rolling_mean) / self->window_count);
  } else {
    self->rolling_mean = reduction_function;
  }

  const float adapted_threshold =
      (TP_UPPER_LIMIT - transient_threshold) * self->rolling_mean;

  memcpy(self->previous_spectrum, spectrum,
         sizeof(float) * (self->half_fft_size + 1));

  if (reduction_function > adapted_threshold) {
    return true;
  }
  return false;
}

static float spectral_flux(const float *spectrum, float *spectrum_prev,
                           const uint32_t half_fft_size) {
  float spectral_flux = 0.f;

  for (uint32_t i = 1; i <= half_fft_size; i++) {
    const float temp = sqrtf(spectrum[i]) - sqrtf(spectrum_prev[i]);
    spectral_flux += (temp + fabsf(temp)) / 2.f;
  }
  return spectral_flux;
}