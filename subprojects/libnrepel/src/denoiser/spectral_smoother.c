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

#include "spectral_smoother.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void get_release_coefficient(SpectralSmoother *self, float release);
static void apply_time_envelope(SpectralSmoother *self);

struct SpectralSmoother {
  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;
  uint32_t hop;

  float *noise_spectrum;

  float *smoothed_spectrum;
  float *smoothed_spectrum_previous;

  float release_coefficient;
};

SpectralSmoother *spectral_smoothing_initialize(const uint32_t fft_size,
                                                const uint32_t sample_rate,
                                                const uint32_t hop) {
  SpectralSmoother *self =
      (SpectralSmoother *)calloc(1U, sizeof(SpectralSmoother));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2U;
  self->sample_rate = sample_rate;
  self->hop = hop;

  self->noise_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));
  self->smoothed_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));
  self->smoothed_spectrum_previous =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));

  self->release_coefficient = 0.f;

  return self;
}

void spectral_smoothing_free(SpectralSmoother *self) {
  free(self->noise_spectrum);
  free(self->smoothed_spectrum);
  free(self->smoothed_spectrum_previous);
  free(self);
}

bool spectral_smoothing_run(SpectralSmoother *self, const float release,
                            const float *signal_spectrum) {
  if (!self || !signal_spectrum) {
    return false;
  }
  get_release_coefficient(self, release);

  memcpy(self->smoothed_spectrum, signal_spectrum,
         sizeof(float) * (self->half_fft_size + 1U));

  apply_time_envelope(self);

  memcpy(self->smoothed_spectrum_previous, self->smoothed_spectrum,
         sizeof(float) * (self->half_fft_size + 1U));

  return true;
}

static void get_release_coefficient(SpectralSmoother *self,
                                    const float release) {
  if (release != 0.f) {
    self->release_coefficient =
        expf(-1000.f / (((release)*self->sample_rate) / self->hop));
  } else {
    self->release_coefficient = 0.f;
  }
}

static void apply_time_envelope(SpectralSmoother *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->smoothed_spectrum[k] > self->smoothed_spectrum_previous[k]) {
      self->smoothed_spectrum[k] =
          self->release_coefficient * self->smoothed_spectrum_previous[k] +
          (1.f - self->release_coefficient) * self->smoothed_spectrum[k];
    }
  }
}