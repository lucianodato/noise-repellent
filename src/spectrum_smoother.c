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

#include "spectrum_smoother.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct SpectralSmoother {
  int fft_size;
  int half_fft_size;
  int samp_rate;
  int hop;

  float *noise_spectrum;
  float *signal_spectrum;

  float *smoothed_spectrum;
  float *smoothed_spectrum_previous;

  float release_coefficient;
};

SpectralSmoother *spectral_smoothing_initialize(int fft_size, int samp_rate,
                                                int hop) {
  SpectralSmoother *self =
      (SpectralSmoother *)calloc(1, sizeof(SpectralSmoother));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->samp_rate = samp_rate;
  self->hop = hop;

  self->signal_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->noise_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->smoothed_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->smoothed_spectrum_previous =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->release_coefficient = 0.f;

  return self;
}

void spectral_smoothing_free(SpectralSmoother *self) {
  free(self->noise_spectrum);
  free(self->signal_spectrum);
  free(self->smoothed_spectrum);
  free(self->smoothed_spectrum_previous);
  free(self);
}

void get_release_coefficient(SpectralSmoother *self, float release) {
  if (release != 0.f) {
    self->release_coefficient =
        expf(-1000.f / (((release)*self->samp_rate) / self->hop));
  } else {
    self->release_coefficient = 0.f;
  }
}

void apply_time_envelope(SpectralSmoother *self) {
  for (int k = 1; k <= self->half_fft_size; k++) {
    if (self->smoothed_spectrum[k] > self->smoothed_spectrum_previous[k]) {
      self->smoothed_spectrum[k] =
          self->release_coefficient * self->smoothed_spectrum_previous[k] +
          (1.f - self->release_coefficient) * self->smoothed_spectrum[k];
    }
  }
}

void spectral_smoothing_run(SpectralSmoother *self, float release) {
  get_release_coefficient(self, release);

  memcpy(self->smoothed_spectrum, self->signal_spectrum,
         sizeof(float) * (self->half_fft_size + 1));

  apply_time_envelope(self);

  memcpy(self->smoothed_spectrum_previous, self->smoothed_spectrum,
         sizeof(float) * (self->half_fft_size + 1));
}