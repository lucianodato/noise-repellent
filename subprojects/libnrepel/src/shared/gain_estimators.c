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

#include "gain_estimators.h"
#include "configurations.h"
#include "general_utils.h"
#include <float.h>
#include <math.h>

void denoise_mixer(const uint32_t fft_size, const uint32_t half_fft_size,
                   float *fft_spectrum, const float *gain_spectrum,
                   float *denoised_spectrum, float *residual_spectrum,
                   const bool residual_listen, const float reduction_amount) {

  // Get denoised spectrum - Apply to both real and complex parts
  for (uint32_t k = 1U; k < half_fft_size; k++) {
    denoised_spectrum[k] = fft_spectrum[k] * gain_spectrum[k];
    denoised_spectrum[fft_size - k] =
        fft_spectrum[fft_size - k] * gain_spectrum[k];
  }

  // Get residual spectrum - Apply to both real and complex parts
  for (uint32_t k = 1U; k < half_fft_size; k++) {
    residual_spectrum[k] = fft_spectrum[k] - denoised_spectrum[k];
    residual_spectrum[fft_size - k] =
        fft_spectrum[fft_size - k] - denoised_spectrum[fft_size - k];
  }

  // Mix denoised and residual
  if (residual_listen) {
    for (uint32_t k = 1U; k < fft_size; k++) {
      fft_spectrum[k] = residual_spectrum[k];
    }
  } else {
    for (uint32_t k = 1U; k < fft_size; k++) {
      fft_spectrum[k] =
          denoised_spectrum[k] + residual_spectrum[k] * reduction_amount;
    }
  }
}

void wiener_subtraction(const uint32_t real_spectrum_size,
                        const float *spectrum, float *gain_spectrum,
                        const float *noise_spectrum) {
  for (uint32_t k = 1U; k < real_spectrum_size; k++) {
    if (noise_spectrum[k] > FLT_MIN) {
      if (spectrum[k] > noise_spectrum[k]) {
        gain_spectrum[k] = (spectrum[k] - noise_spectrum[k]) / spectrum[k];
      } else {
        gain_spectrum[k] = 0.F;
      }
    } else {
      gain_spectrum[k] = 1.F;
    }
  }
}

void spectral_gating(const uint32_t real_spectrum_size, const float *spectrum,
                     float *gain_spectrum, const float *noise_spectrum) {
  for (uint32_t k = 1U; k < real_spectrum_size; k++) {
    if (noise_spectrum[k] > FLT_MIN) {
      if (spectrum[k] >= noise_spectrum[k]) {
        gain_spectrum[k] = 1.F;
      } else {
        gain_spectrum[k] = 0.F;
      }
    } else {
      gain_spectrum[k] = 1.F;
    }
  }
}

void generalized_spectral_subtraction(const uint32_t real_spectrum_size,
                                      const float *alpha, const float *beta,
                                      const float *spectrum,
                                      const float *noise_spectrum,
                                      float *gain_spectrum) {
  for (uint32_t k = 1U; k < real_spectrum_size; k++) {
    if (spectrum[k] > FLT_MIN) {
      if (powf((noise_spectrum[k] / spectrum[k]), GAMMA1) <
          (1.F / (alpha[k] + beta[k]))) {
        gain_spectrum[k] =
            fmaxf(powf(1.F - (alpha[k] *
                              powf((noise_spectrum[k] / spectrum[k]), GAMMA1)),
                       GAMMA2),
                  0.F);
      } else {
        gain_spectrum[k] = fmaxf(
            powf(beta[k] * powf((noise_spectrum[k] / spectrum[k]), GAMMA1),
                 GAMMA2),
            0.F);
      }
    } else {
      gain_spectrum[k] = 1.F;
    }
  }
}