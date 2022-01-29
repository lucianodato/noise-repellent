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

#include "spectral_denoiser.h"
#include "../shared/configurations.h"
#include "../shared/gain_estimators.h"
#include "../shared/noise_estimator.h"
#include "../shared/oversubtraction_criterias.h"
#include "../shared/postfilter.h"
#include "../shared/spectral_features.h"
#include "../shared/spectral_smoother.h"
#include "../shared/spectral_utils.h"
#include "../shared/spectral_whitening.h"
#include "../shared/transient_detector.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct SbSpectralDenoiser {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;
  float default_oversubtraction;
  float default_undersubtraction;
  bool transient_detected;

  float *gain_spectrum;
  float *residual_spectrum;
  float *denoised_spectrum;

  SpectrumType spectrum_type;
  OversubtractionType oversubtraction_type;
  CriticalBandType band_type;
  DenoiserParameters denoise_parameters;

  NoiseEstimator *noise_estimator;
  SpectralWhitening *whitener;
  PostFilter *postfiltering;
  NoiseProfile *noise_profile;
  SpectralFeatures *spectral_features;

  OversubtractionCriterias *oversubtraction_criteria;
  TransientDetector *transient_detection;
  SpectralSmoother *spectrum_smoothing;
} SbSpectralDenoiser;

SpectralProcessorHandle spectral_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile *noise_profile) {

  SbSpectralDenoiser *self =
      (SbSpectralDenoiser *)calloc(1U, sizeof(SbSpectralDenoiser));

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->hop = self->fft_size / overlap_factor;
  self->sample_rate = sample_rate;
  self->spectrum_type = SPECTRAL_TYPE_GENERAL;
  self->oversubtraction_type = OVERSUBTRACTION_TYPE;
  self->band_type = CRITICAL_BANDS_TYPE;
  self->default_oversubtraction = DEFAULT_OVERSUBTRACTION;
  self->default_undersubtraction = DEFAULT_UNDERSUBTRACTION;

  self->gain_spectrum = (float *)calloc(self->fft_size, sizeof(float));
  initialize_spectrum_with_value(self->gain_spectrum, self->fft_size, 1.F);

  self->noise_profile = noise_profile;

  self->noise_estimator =
      noise_estimation_initialize(self->fft_size, noise_profile);

  self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
  self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));

  self->spectral_features =
      spectral_features_initialize(self->real_spectrum_size);

  self->postfiltering = postfilter_initialize(self->fft_size);

  self->transient_detection = transient_detector_initialize(self->fft_size);
  self->spectrum_smoothing = spectral_smoothing_initialize(
      self->fft_size, self->sample_rate, self->hop);

  self->oversubtraction_criteria = oversubtraction_criterias_initialize(
      self->oversubtraction_type, self->fft_size, self->band_type,
      self->sample_rate, self->spectrum_type);

  self->whitener = spectral_whitening_initialize(self->fft_size,
                                                 self->sample_rate, self->hop);

  return self;
}

void spectral_denoiser_free(SpectralProcessorHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  spectral_whitening_free(self->whitener);
  noise_estimation_free(self->noise_estimator);
  spectral_features_free(self->spectral_features);
  transient_detector_free(self->transient_detection);
  spectral_smoothing_free(self->spectrum_smoothing);
  oversubtraction_criterias_free(self->oversubtraction_criteria);
  postfilter_free(self->postfiltering);

  free(self->residual_spectrum);
  free(self->denoised_spectrum);
  free(self->gain_spectrum);
  free(self);
}

bool load_reduction_parameters(SpectralProcessorHandle instance,
                               DenoiserParameters parameters) {
  if (!instance) {
    return false;
  }

  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;
  self->denoise_parameters = parameters;

  return true;
}

bool spectral_denoiser_run(SpectralProcessorHandle instance,
                           float *fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  float *reference_spectrum =
      get_spectral_feature(self->spectral_features, fft_spectrum,
                           self->fft_size, self->spectrum_type);

  if (self->denoise_parameters.learn_noise) {
    noise_estimation_run(self->noise_estimator, reference_spectrum);
  }

  if (is_noise_estimation_available(self->noise_profile)) {
    float *noise_profile = get_noise_profile(self->noise_profile);

    if (self->denoise_parameters.transient_threshold > 1.F) {
      self->transient_detected = transient_detector_run(
          self->transient_detection,
          self->denoise_parameters.transient_threshold, reference_spectrum);
    }

    OversustractionParameters oversubtraction_parameters =
        (OversustractionParameters){
            .oversubtraction = self->default_oversubtraction +
                               self->denoise_parameters.noise_rescale,
            .undersubtraction = self->default_undersubtraction,
        };
    apply_oversustraction_criteria(self->oversubtraction_criteria,
                                   reference_spectrum, noise_profile,
                                   oversubtraction_parameters);

    spectral_smoothing_run(self->spectrum_smoothing,
                           self->denoise_parameters.release_time,
                           reference_spectrum);

    if (self->transient_detected &&
        self->denoise_parameters.transient_threshold > 1.F) {
      wiener_subtraction(self->real_spectrum_size, self->fft_size,
                         reference_spectrum, self->gain_spectrum,
                         noise_profile);
    } else {
      spectral_gating(self->real_spectrum_size, self->fft_size,
                      reference_spectrum, self->gain_spectrum, noise_profile);
    }

    // Apply post filtering to reduce residual noise on low SNR frames
    postfilter_apply(self->postfiltering, fft_spectrum, self->gain_spectrum);

    // FIXME (luciano/fix): Apply whitening to gain weights instead of the
    // resulting noise residue
    if (self->denoise_parameters.whitening_factor > 0.F) {
      spectral_whitening_run(self->whitener,
                             self->denoise_parameters.whitening_factor,
                             self->gain_spectrum);
    }

    denoise_mixer(self->fft_size, fft_spectrum, self->gain_spectrum,
                  self->denoised_spectrum, self->residual_spectrum,
                  self->denoise_parameters.residual_listen,
                  self->denoise_parameters.reduction_amount);
  }

  return true;
}