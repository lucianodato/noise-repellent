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

#include "noise_estimator.h"
#include "../../shared/configurations.h"
#include "../../shared/spectral_features.h"
#include "louizou_estimator.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseEstimator {
  uint32_t fft_size;
  uint32_t half_fft_size;
  bool noise_spectrum_available;

  NoiseProfile *noise_profile;
  LouizouEstimator *adaptive_estimator;
};

NoiseEstimator *noise_estimation_initialize(const uint32_t fft_size,
                                            const uint32_t sample_rate,
                                            NoiseProfile *noise_profile) {
  NoiseEstimator *self = (NoiseEstimator *)calloc(1U, sizeof(NoiseEstimator));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2U;

  self->noise_spectrum_available = false;

  self->noise_profile = noise_profile;

  self->adaptive_estimator = louizou_estimator_initialize(
      self->half_fft_size + 1U, sample_rate, fft_size);

  return self;
}

void noise_estimation_free(NoiseEstimator *self) {
  louizou_estimator_free(self->adaptive_estimator);
  free(self);
}

bool is_noise_estimation_available(NoiseEstimator *self) {
  return self->noise_spectrum_available;
}

bool noise_estimation_run(NoiseEstimator *self, float *signal_spectrum) {
  if (!self || !signal_spectrum) {
    return false;
  }

  increment_blocks_averaged(self->noise_profile);

  float *noise_profile = get_noise_profile(self->noise_profile);

  get_rolling_mean_spectrum(
      noise_profile, signal_spectrum,
      get_noise_profile_blocks_averaged(self->noise_profile),
      self->half_fft_size);

  if (get_noise_profile_blocks_averaged(self->noise_profile) >
      MIN_NUMBER_OF_WINDOWS_NOISE_AVERAGED) {
    self->noise_spectrum_available = true;
  }

  return true;
}

bool noise_estimation_run_adaptive(NoiseEstimator *self,
                                   float *signal_spectrum) {
  if (!self || !signal_spectrum) {
    return false;
  }

  float *noise_profile = get_noise_profile(self->noise_profile);

  louizou_estimator_run(self->adaptive_estimator, signal_spectrum,
                        noise_profile);

  self->noise_spectrum_available = true;

  return true;
}
