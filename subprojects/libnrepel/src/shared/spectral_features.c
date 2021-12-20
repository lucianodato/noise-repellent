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

#include "spectral_features.h"
#include <math.h>
#include <stdlib.h>

struct SpectralFeatures {
  float *power_spectrum;
  float *phase_spectrum;
  float *magnitude_spectrum;

  uint32_t spectral_size;
};

SpectralFeatures *spectral_features_initialize(const uint32_t spectral_size) {
  SpectralFeatures *self =
      (SpectralFeatures *)calloc(1, sizeof(SpectralFeatures));

  self->spectral_size = spectral_size;

  self->power_spectrum =
      (float *)calloc((self->spectral_size + 1), sizeof(float));
  self->phase_spectrum =
      (float *)calloc((self->spectral_size + 1), sizeof(float));
  self->magnitude_spectrum =
      (float *)calloc((self->spectral_size + 1), sizeof(float));

  return self;
}

void spectral_features_free(SpectralFeatures *self) {
  free(self->power_spectrum);
  free(self->phase_spectrum);
  free(self->magnitude_spectrum);
  free(self);
}

float *get_power_spectrum(SpectralFeatures *self) {
  return self->power_spectrum;
}
float *get_magnitude_spectrum(SpectralFeatures *self) {
  return self->magnitude_spectrum;
}
float *get_phase_spectrum(SpectralFeatures *self) {
  return self->phase_spectrum;
}

bool compute_power_spectrum(SpectralFeatures *self, const float *fft_spectrum,
                            const uint32_t fft_spectrum_size) {
  if (!self || !fft_spectrum || !fft_spectrum_size) {
    return false;
  }

  float real_bin = fft_spectrum[0];

  self->power_spectrum[0] = real_bin * real_bin;

  for (uint32_t k = 1; k <= self->spectral_size; k++) {
    float power = 0.f;

    real_bin = fft_spectrum[k];
    float imag_bin = fft_spectrum[fft_spectrum_size - k];

    if (k < self->spectral_size) {
      power = (real_bin * real_bin + imag_bin * imag_bin);
    } else {
      power = real_bin * real_bin;
    }

    self->power_spectrum[k] = power;
  }

  return true;
}

bool compute_magnitude_spectrum(SpectralFeatures *self,
                                const float *fft_spectrum,
                                const uint32_t fft_spectrum_size) {
  if (!self || !fft_spectrum || !fft_spectrum_size) {
    return false;
  }

  float real_bin = fft_spectrum[0];

  self->magnitude_spectrum[0] = real_bin;

  for (uint32_t k = 1; k <= self->spectral_size; k++) {
    float magnitude = 0.f;

    real_bin = fft_spectrum[k];
    float imag_bin = fft_spectrum[fft_spectrum_size - k];

    if (k < self->spectral_size) {
      magnitude = sqrtf(real_bin * real_bin + imag_bin * imag_bin);

    } else {
      magnitude = real_bin;
    }

    self->magnitude_spectrum[k] = magnitude;
  }

  return true;
}

bool compute_phase_spectrum(SpectralFeatures *self, const float *fft_spectrum,
                            const uint32_t fft_spectrum_size) {
  if (!self || !fft_spectrum || !fft_spectrum_size) {
    return false;
  }

  float real_bin = fft_spectrum[0];
  self->phase_spectrum[0] = atan2f(real_bin, 0.f);

  for (uint32_t k = 1; k <= self->spectral_size; k++) {
    float phase = 0.f;

    real_bin = fft_spectrum[k];
    float imag_bin = fft_spectrum[fft_spectrum_size - k];

    if (k < self->spectral_size) {
      phase = atan2f(real_bin, imag_bin);
    } else {
      phase = atan2f(real_bin, 0.f);
    }

    self->phase_spectrum[k] = phase;
  }

  return true;
}