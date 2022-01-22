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
#include "../shared/adaptive_noise_estimator.h"
#include "../shared/configurations.h"
#include "../shared/gain_estimators.h"
#include "../shared/oversubtraction_criterias.h"
#include "../shared/spectral_features.h"
#include "../shared/spectral_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct SpectralAdaptiveDenoiser {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  float default_masking_ceiling;
  float default_masking_floor;

  AdaptiveDenoiserParameters parameters;

  float *gain_spectrum;
  float *residual_spectrum;
  float *denoised_spectrum;
  float *noise_profile;

  SpectrumType spectrum_type;

  OversubtractionCriterias *oversubtraction_criteria;
  AdaptiveNoiseEstimator *adaptive_estimator;
  SpectralFeatures *spectral_features;
} SpectralAdaptiveDenoiser;

SpectralProcessorHandle
spectral_adaptive_denoiser_initialize(const uint32_t sample_rate,
                                      const uint32_t fft_size) {

  SpectralAdaptiveDenoiser *self =
      (SpectralAdaptiveDenoiser *)calloc(1U, sizeof(SpectralAdaptiveDenoiser));

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->sample_rate = sample_rate;
  self->default_masking_ceiling = DEFAULT_MASKING_CEILING;
  self->default_masking_floor = DEFAULT_MASKING_FLOOR;
  self->spectrum_type = SPECTRAL_TYPE_SPEECH;

  self->gain_spectrum =
      (float *)calloc(self->real_spectrum_size, sizeof(float));
  initialize_spectrum_with_value(self->gain_spectrum, self->real_spectrum_size,
                                 1.F);
  self->noise_profile =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

  self->adaptive_estimator = louizou_estimator_initialize(
      self->real_spectrum_size, sample_rate, fft_size);

  self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
  self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));

  self->oversubtraction_criteria = oversubtraction_criterias_initialize(
      MASKING_THRESHOLDS, N_CRITICAL_BANDS_SPEECH, self->fft_size,
      CRITICAL_BANDS_TYPE_SPEECH, self->sample_rate, self->spectrum_type);

  self->spectral_features =
      spectral_features_initialize(self->real_spectrum_size);

  return self;
}

void spectral_adaptive_denoiser_free(SpectralProcessorHandle instance) {
  SpectralAdaptiveDenoiser *self = (SpectralAdaptiveDenoiser *)instance;

  louizou_estimator_free(self->adaptive_estimator);
  spectral_features_free(self->spectral_features);
  oversubtraction_criterias_free(self->oversubtraction_criteria);

  free(self->residual_spectrum);
  free(self->denoised_spectrum);
  free(self->gain_spectrum);
  free(self);
}

bool load_adaptive_reduction_parameters(SpectralProcessorHandle instance,
                                        AdaptiveDenoiserParameters parameters) {
  if (!instance) {
    return false;
  }

  SpectralAdaptiveDenoiser *self = (SpectralAdaptiveDenoiser *)instance;
  self->parameters = parameters;

  return true;
}

bool spectral_adaptive_denoiser_run(SpectralProcessorHandle instance,
                                    float *fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SpectralAdaptiveDenoiser *self = (SpectralAdaptiveDenoiser *)instance;

  float *reference_spectrum =
      get_spectral_feature(self->spectral_features, fft_spectrum,
                           self->fft_size, self->spectrum_type);

  // Estimate noise
  louizou_estimator_run(self->adaptive_estimator, reference_spectrum,
                        self->noise_profile);

  // Scale estimated noise profile for oversubtraction using a posteriori SNR
  // piecewise function
  OversustractionParameters oversubtraction_parameters =
      (OversustractionParameters){
          .noise_rescale = self->parameters.noise_rescale,
          .masking_ceil = self->default_masking_ceiling,
          .masking_floor = self->default_masking_floor,
      };
  apply_oversustraction_criteria(self->oversubtraction_criteria,
                                 reference_spectrum, self->noise_profile,
                                 oversubtraction_parameters);

  // Get reduction gain weights
  wiener_subtraction(self->real_spectrum_size, reference_spectrum,
                     self->gain_spectrum, self->noise_profile);

  // Mix results
  denoise_mixer(self->fft_size, self->real_spectrum_size, fft_spectrum,
                self->gain_spectrum, self->denoised_spectrum,
                self->residual_spectrum, self->parameters.residual_listen,
                self->parameters.reduction_amount);

  return true;
}