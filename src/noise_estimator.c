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
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseEstimator {
  int fft_size;
  int half_fft_size;
  bool noise_spectrum_available;
  float noise_blocks_count;
};

NoiseEstimator *noise_estimation_initialize(const int fft_size) {
  NoiseEstimator *self = (NoiseEstimator *)calloc(1, sizeof(NoiseEstimator));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->noise_blocks_count = 0;
  self->noise_spectrum_available = false;

  return self;
}

void noise_estimation_free(NoiseEstimator *self) { free(self); }

bool is_noise_estimation_available(NoiseEstimator *self) {
  return self->noise_spectrum_available;
}

void noise_estimation_run(NoiseEstimator *self, float *noise_spectrum,
                          const float *spectrum) {
  self->noise_blocks_count++;

  for (int k = 1; k <= self->half_fft_size; k++) {
    if (self->noise_blocks_count <= 1.f) {
      noise_spectrum[k] = spectrum[k];
    } else {
      noise_spectrum[k] +=
          ((spectrum[k] - noise_spectrum[k]) / self->noise_blocks_count);
    }
  }

  self->noise_spectrum_available = true;
}
