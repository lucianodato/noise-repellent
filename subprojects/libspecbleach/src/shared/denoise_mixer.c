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

#include "denoise_mixer.h"
#include "../shared/spectral_whitening.h"
#include <stdlib.h>
#include <string.h>

struct DenoiseMixer {
  SpectralWhitening *whitener;

  float *residual_spectrum;
  float *denoised_spectrum;

  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;
};

DenoiseMixer *denoise_mixer_initialize(uint32_t fft_size, uint32_t sample_rate,
                                       uint32_t hop) {
  DenoiseMixer *self = (DenoiseMixer *)calloc(1U, sizeof(DenoiseMixer));

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->sample_rate = sample_rate;
  self->hop = hop;

  self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
  self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));

  self->whitener = spectral_whitening_initialize(self->fft_size,
                                                 self->sample_rate, self->hop);

  return self;
}

void denoise_mixer_free(DenoiseMixer *self) {
  spectral_whitening_free(self->whitener);
  free(self->residual_spectrum);
  free(self->denoised_spectrum);
  free(self);
}

bool denoise_mixer_run(DenoiseMixer *self, float *fft_spectrum,
                       const float *gain_spectrum,
                       DenoiseMixerParameters parameters) {

  if (!fft_spectrum || !gain_spectrum) {
    return false;
  }

  // Get denoised spectrum - Apply to both real and complex parts
  for (uint32_t k = 1U; k < self->fft_size; k++) {
    self->denoised_spectrum[k] = fft_spectrum[k] * gain_spectrum[k];
  }

  // Get residual spectrum - Apply to both real and complex parts
  for (uint32_t k = 1U; k < self->fft_size; k++) {
    self->residual_spectrum[k] = fft_spectrum[k] - self->denoised_spectrum[k];
  }

  if (parameters.whitening_amount > 0.F) {
    spectral_whitening_run(self->whitener, parameters.whitening_amount,
                           self->residual_spectrum);
  }

  // Mix denoised and residual
  if (parameters.residual_listen) {
    for (uint32_t k = 1U; k < self->fft_size; k++) {
      fft_spectrum[k] = self->residual_spectrum[k];
    }
  } else {
    for (uint32_t k = 1U; k < self->fft_size; k++) {
      fft_spectrum[k] = self->denoised_spectrum[k] +
                        self->residual_spectrum[k] * parameters.noise_level;
    }
  }

  return true;
}