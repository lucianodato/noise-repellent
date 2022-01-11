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

#include "gain_estimator.h"
#include "../../shared/configurations.h"
#include "../../shared/spectral_utils.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct GainEstimatorAdaptive {
  uint32_t fft_size;
  uint32_t half_fft_size;

  NrepelDenoiseParameters *denoise_parameters;
  NoiseProfile *noise_profile;
};

GainEstimatorAdaptive *
adaptive_gain_estimation_initialize(const uint32_t fft_size,
                                    NrepelDenoiseParameters *parameters,
                                    NoiseProfile *noise_profile) {
  GainEstimatorAdaptive *self =
      (GainEstimatorAdaptive *)calloc(1U, sizeof(GainEstimatorAdaptive));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2U;

  self->noise_profile = noise_profile;
  self->denoise_parameters = parameters;

  return self;
}

void adaptive_gain_estimation_free(GainEstimatorAdaptive *self) { free(self); }

bool gain_estimation_run_adaptive(GainEstimatorAdaptive *self,
                                  const float *signal_spectrum,
                                  float *gain_spectrum) {
  if (!self || !signal_spectrum || !gain_spectrum) {
    return false;
  }

  float *noise_profile = get_noise_profile(self->noise_profile);

  for (uint32_t k = 1U; k <= self->half_fft_size; k++) {
    noise_profile[k] =
        noise_profile[k] * self->denoise_parameters->noise_rescale;
  }

  wiener_subtraction(self->half_fft_size, signal_spectrum, gain_spectrum,
                     noise_profile);

  return true;
}
