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

#include "noise_estimator.h"
#include "../shared/common.h"
#include "../shared/spectral_features.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseEstimator {
  uint32_t fft_size;
  uint32_t half_fft_size;
  bool noise_spectrum_available;
  float noise_blocks_count;

  SpectralFeatures processing_spectrums;

  NoiseProfile *noise_profile;
};

NoiseEstimator *noise_estimation_initialize(const uint32_t fft_size,
                                            NoiseProfile *noise_profile) {
  NoiseEstimator *self = (NoiseEstimator *)calloc(1, sizeof(NoiseEstimator));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->noise_blocks_count = 0;
  self->noise_spectrum_available = false;

  self->noise_profile = noise_profile;

  self->processing_spectrums.power_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  return self;
}

void noise_estimation_free(NoiseEstimator *self) {
  free(self->processing_spectrums.power_spectrum);
  free(self);
}

bool is_noise_estimation_available(NoiseEstimator *self) {
  return self->noise_spectrum_available;
}

void noise_estimation_run(SPECTRAL_PROCESSOR instance, float *fft_spectrum) {
  NoiseEstimator *self = (NoiseEstimator *)instance;

  self->noise_blocks_count++;

  get_fft_power_spectrum(fft_spectrum, self->fft_size,
                         self->processing_spectrums.power_spectrum,
                         self->half_fft_size);

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->noise_blocks_count <= 1.f) {
      get_noise_profile(self->noise_profile)[k] =
          self->processing_spectrums.power_spectrum[k];
    } else {
      get_noise_profile(self->noise_profile)[k] +=
          ((self->processing_spectrums.power_spectrum[k] -
            get_noise_profile(self->noise_profile)[k]) /
           self->noise_blocks_count);
    }
  }

  self->noise_spectrum_available = true;
}
