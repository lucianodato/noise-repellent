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

struct TransientDetector {
  int fft_size;
  int half_fft_size;

  float *spectrum;

  float *previous_spectrum;
  float r_mean;
  bool transient_present;
  float window_count;
};

TransientDetector *transient_detector_initialize(int fft_size) {
  TransientDetector *self =
      (TransientDetector *)calloc(1, sizeof(TransientDetector));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;

  self->spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->previous_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->window_count = 0.f;
  self->r_mean = 0.f;
  self->transient_present = false;

  return self;
}

void transient_detector_free(TransientDetector *self) {
  free(self->spectrum);
  free(self->previous_spectrum);
  free(self);
}

float spectral_flux(float *spectrum, float *spectrum_prev, float N) {
  int i = 0;
  float spectral_flux = 0.f;

  for (i = 1; i <= N; i++) {
    float temp = 0.f;
    temp = sqrtf(spectrum[i]) - sqrtf(spectrum_prev[i]);
    spectral_flux += (temp + fabsf(temp)) / 2.f;
  }
  return spectral_flux;
}

bool transient_detector_run(TransientDetector *self,
                            float transient_threshold) {
  float adapted_threshold = 0.f;
  float reduction_function = 0.f;

  reduction_function = spectral_flux(self->spectrum, self->previous_spectrum,
                                     self->half_fft_size);

  self->window_count += 1.f;

  if (self->window_count > 1.f) {
    self->r_mean += ((reduction_function - self->r_mean) / self->window_count);
  } else {
    self->r_mean = reduction_function;
  }

  adapted_threshold = (TP_UPPER_LIMIT - transient_threshold) * self->r_mean;

  memcpy(self->previous_spectrum, self->spectrum,
         sizeof(float) * (self->half_fft_size + 1));

  if (reduction_function > adapted_threshold) {
    return true;
  }
  return false;
}