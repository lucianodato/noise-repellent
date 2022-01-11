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

#include "adaptive_denoiser.h"
#include "../shared/configurations.h"
#include "../shared/spectral_features.h"
#include "../shared/spectral_utils.h"
#include "gainmodel/gain_estimator.h"
#include "noisemodel/noise_estimator.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct SpectralAdaptiveDenoiser {
  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;

  float *gain_spectrum;
  float *residual_spectrum;
  float *denoised_spectrum;

  NoiseEstimatorAdaptive *noise_estimator;
  GainEstimatorAdaptive *gain_estimator;
  NrepelDenoiseParameters *denoise_parameters;
  SpectralFeatures *spectral_features;
} SpectralAdaptiveDenoiser;

SpectralProcessorHandle spectral_adaptive_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    NoiseProfile *noise_profile, NrepelDenoiseParameters *parameters) {

  SpectralAdaptiveDenoiser *self =
      (SpectralAdaptiveDenoiser *)calloc(1U, sizeof(SpectralAdaptiveDenoiser));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2U;
  self->sample_rate = sample_rate;

  self->gain_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));
  initialize_spectrum_with_value(self->gain_spectrum, self->half_fft_size + 1U,
                                 1.F);

  self->denoise_parameters = parameters;

  self->noise_estimator = adaptive_noise_estimation_initialize(
      self->fft_size, sample_rate, noise_profile);

  self->gain_estimator = adaptive_gain_estimation_initialize(
      self->fft_size, self->denoise_parameters, noise_profile);

  self->residual_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));
  self->denoised_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));

  self->spectral_features =
      spectral_features_initialize(self->half_fft_size + 1U);

  return self;
}

void spectral_adaptive_denoiser_free(SpectralProcessorHandle instance) {
  SpectralAdaptiveDenoiser *self = (SpectralAdaptiveDenoiser *)instance;

  adaptive_gain_estimation_free(self->gain_estimator);
  adaptive_noise_estimation_free(self->noise_estimator);
  spectral_features_free(self->spectral_features);

  free(self->residual_spectrum);
  free(self->denoised_spectrum);
  free(self->gain_spectrum);
  free(self);
}

bool spectral_adaptive_denoiser_run(SpectralProcessorHandle instance,
                                    float *fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SpectralAdaptiveDenoiser *self = (SpectralAdaptiveDenoiser *)instance;

  float *reference_spectrum = get_spectral_feature(
      self->spectral_features, fft_spectrum, self->fft_size, SPECTRAL_TYPE);

  noise_estimation_run_adaptive(self->noise_estimator, reference_spectrum);
  gain_estimation_run_adaptive(self->gain_estimator, reference_spectrum,
                               self->gain_spectrum);
  denoise_mixer(self->half_fft_size, fft_spectrum, self->gain_spectrum,
                self->denoised_spectrum, self->residual_spectrum,
                self->denoise_parameters->residual_listen,
                self->denoise_parameters->reduction_amount);

  return true;
}