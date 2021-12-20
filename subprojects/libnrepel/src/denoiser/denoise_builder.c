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

#include "denoise_builder.h"
#include "spectral_whitening.h"
#include <stdlib.h>
#include <string.h>

struct DenoiseBuilder {
  float *residual_spectrum;
  float *denoised_spectrum;

  uint32_t half_fft_size;
  SpectralWhitening *whitener;
};

DenoiseBuilder *denoise_builder_initialize(const uint32_t half_fft_size,
                                           const uint32_t fft_size,
                                           const uint32_t hop,
                                           const uint32_t sample_rate) {
  DenoiseBuilder *self = (DenoiseBuilder *)calloc(1, sizeof(DenoiseBuilder));

  self->half_fft_size = half_fft_size;

  self->residual_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->denoised_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->whitener = spectral_whitening_initialize(fft_size, sample_rate, hop);

  return self;
}

void denoise_builder_free(DenoiseBuilder *self) {
  spectral_whitening_free(self->whitener);

  free(self->residual_spectrum);
  free(self->denoised_spectrum);
  free(self);
}

static void get_denoised_spectrum(DenoiseBuilder *self,
                                  const float *fft_spectrum,
                                  const float *gain_spectrum) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->denoised_spectrum[k] = fft_spectrum[k] * gain_spectrum[k];
  }
}

static void get_residual_spectrum(DenoiseBuilder *self,
                                  const float *fft_spectrum,
                                  const float whitening_factor) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->residual_spectrum[k] = fft_spectrum[k] - self->denoised_spectrum[k];
  }

  if (whitening_factor > 0.f) {
    spectral_whitening_run(self->whitener, whitening_factor,
                           self->residual_spectrum);
  }
}

void denoise_build(DenoiseBuilder *self, const bool residual_listening,
                   const float reduction_amout, const float whitening_factor,
                   const float *gain_spectrum, float *fft_spectrum) {

  get_denoised_spectrum(self, fft_spectrum, gain_spectrum);

  get_residual_spectrum(self, fft_spectrum, whitening_factor);

  if (residual_listening) {
    for (uint32_t k = 1; k <= self->half_fft_size; k++) {
      fft_spectrum[k] = self->residual_spectrum[k];
    }
  } else {
    for (uint32_t k = 1; k <= self->half_fft_size; k++) {
      fft_spectrum[k] = self->denoised_spectrum[k] +
                        self->residual_spectrum[k] * reduction_amout;
    }
  }
}