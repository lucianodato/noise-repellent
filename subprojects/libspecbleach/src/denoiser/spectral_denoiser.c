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
#include "../shared/denoise_mixer.h"
#include "../shared/gain_estimators.h"
#include "../shared/noise_estimator.h"
#include "../shared/noise_scaling_criterias.h"
#include "../shared/postfilter.h"
#include "../shared/spectral_features.h"
#include "../shared/spectral_smoother.h"
#include "../shared/spectral_utils.h"
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

  float *gain_spectrum;
  float *alpha;
  float *beta;
  float *noise_spectrum;

  SpectrumType spectrum_type;
  NoiseScalingType noise_scaling_type;
  CriticalBandType band_type;
  DenoiserParameters denoise_parameters;
  GainEstimationType gain_estimation_type;
  TimeSmoothingType time_smoothing_type;

  NoiseEstimator *noise_estimator;
  PostFilter *postfiltering;
  NoiseProfile *noise_profile;
  SpectralFeatures *spectral_features;
  DenoiseMixer *mixer;

  NoiseScalingCriterias *noise_scaling_criteria;
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
  self->noise_scaling_type = OVERSUBTRACTION_TYPE;
  self->band_type = CRITICAL_BANDS_TYPE;
  self->default_oversubtraction = DEFAULT_OVERSUBTRACTION;
  self->default_undersubtraction = DEFAULT_UNDERSUBTRACTION;
  self->gain_estimation_type = GAIN_ESTIMATION_TYPE;
  self->time_smoothing_type = TIME_SMOOTHING_TYPE;

  self->gain_spectrum = (float *)calloc(self->fft_size, sizeof(float));
  initialize_spectrum_with_value(self->gain_spectrum, self->fft_size, 1.F);
  self->alpha = (float *)calloc(self->real_spectrum_size, sizeof(float));
  initialize_spectrum_with_value(self->alpha, self->real_spectrum_size, 1.F);
  self->beta = (float *)calloc(self->real_spectrum_size, sizeof(float));

  self->noise_profile = noise_profile;
  self->noise_spectrum =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

  self->noise_estimator =
      noise_estimation_initialize(self->fft_size, noise_profile);

  self->spectral_features =
      spectral_features_initialize(self->real_spectrum_size);

  self->postfiltering = postfilter_initialize(self->fft_size);

  self->spectrum_smoothing = spectral_smoothing_initialize(
      self->fft_size, self->sample_rate, self->hop, self->time_smoothing_type);

  self->noise_scaling_criteria = noise_scaling_criterias_initialize(
      self->noise_scaling_type, self->fft_size, self->band_type,
      self->sample_rate, self->spectrum_type);

  self->mixer =
      denoise_mixer_initialize(self->fft_size, self->sample_rate, self->hop);

  return self;
}

void spectral_denoiser_free(SpectralProcessorHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  noise_estimation_free(self->noise_estimator);
  spectral_features_free(self->spectral_features);
  spectral_smoothing_free(self->spectrum_smoothing);
  noise_scaling_criterias_free(self->noise_scaling_criteria);
  postfilter_free(self->postfiltering);
  denoise_mixer_free(self->mixer);

  free(self->gain_spectrum);
  free(self->alpha);
  free(self->beta);
  free(self->noise_spectrum);
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
  } else if (is_noise_estimation_available(self->noise_profile)) {
    memcpy(self->noise_spectrum, get_noise_profile(self->noise_profile),
           self->real_spectrum_size * sizeof(float));

    NoiseScalingParameters oversubtraction_parameters =
        (NoiseScalingParameters){
            .oversubtraction = self->default_oversubtraction +
                               self->denoise_parameters.noise_rescale,
            .undersubtraction = self->default_undersubtraction,
        };
    apply_noise_scaling_criteria(
        self->noise_scaling_criteria, reference_spectrum, self->noise_spectrum,
        self->alpha, self->beta, oversubtraction_parameters);

    spectral_smoothing_run(self->spectrum_smoothing,
                           self->denoise_parameters.release_time,
                           reference_spectrum, self->noise_spectrum);

    // Get reduction gain weights
    estimate_gains(self->real_spectrum_size, self->fft_size, reference_spectrum,
                   self->noise_spectrum, self->gain_spectrum, self->alpha,
                   self->beta, self->gain_estimation_type);

    // Apply post filtering to reduce residual noise on low SNR frames
    postfilter_apply(self->postfiltering, fft_spectrum, self->gain_spectrum);

    DenoiseMixerParameters mixer_parameters = (DenoiseMixerParameters){
        .noise_level = self->denoise_parameters.reduction_amount,
        .residual_listen = self->denoise_parameters.residual_listen,
        .whitening_amount = self->denoise_parameters.whitening_factor,
    };

    denoise_mixer_run(self->mixer, fft_spectrum, self->gain_spectrum,
                      mixer_parameters);
  }

  return true;
}