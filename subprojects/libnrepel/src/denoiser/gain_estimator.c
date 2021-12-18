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

#include "gain_estimator.h"
#include "masking_estimator.h"
#include "spectrum_smoother.h"
#include "transient_detector.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define GAMMA1 2.f
#define GAMMA2 0.5f

#define ALPHA_MAX 6.f
#define ALPHA_MIN 1.f
#define BETA_MAX 0.02f
#define BETA_MIN 0.0f

static void wiener_subtraction(GainEstimator *self);
static void spectral_gating(GainEstimator *self);
static void compute_alpha_and_beta(GainEstimator *self,
                                   float masking_ceiling_limit,
                                   float masking_floor_limit);

struct GainEstimator {
  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;
  uint32_t hop;

  float *gain_spectrum;
  float *noise_profile;
  float *signal_spectrum;
  float *alpha;
  float *beta;
  float *masking_thresholds;
  float *clean_signal_estimation;

  bool transient_detected;

  ProcessorParameters *denoise_parameters;
  MaskingEstimator *masking_estimation;
  TransientDetector *transient_detection;
  SpectralSmoother *spectrum_smoothing;
};

GainEstimator *gain_estimation_initialize(const uint32_t fft_size,
                                          const uint32_t sample_rate,
                                          const uint32_t hop,
                                          ProcessorParameters *parameters) {
  GainEstimator *self = (GainEstimator *)calloc(1, sizeof(GainEstimator));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->sample_rate = sample_rate;
  self->hop = hop;

  self->signal_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->gain_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->noise_profile =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->alpha = (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->beta = (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->masking_thresholds =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->clean_signal_estimation =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->masking_estimation =
      masking_estimation_initialize(self->fft_size, self->sample_rate);
  self->transient_detection = transient_detector_initialize(self->fft_size);
  self->spectrum_smoothing = spectral_smoothing_initialize(
      self->fft_size, self->sample_rate, self->hop);

  self->denoise_parameters = parameters;

  return self;
}

void gain_estimation_free(GainEstimator *self) {
  masking_estimation_free(self->masking_estimation);
  transient_detector_free(self->transient_detection);
  spectral_smoothing_free(self->spectrum_smoothing);

  free(self->noise_profile);
  free(self->gain_spectrum);
  free(self->signal_spectrum);
  free(self->alpha);
  free(self->beta);
  free(self->masking_thresholds);
  free(self->clean_signal_estimation);
  free(self);
}

void gain_estimation_run(GainEstimator *self, const float *signal_spectrum,
                         const float *noise_profile, float *gain_spectrum) {
  memcpy(self->signal_spectrum, signal_spectrum,
         sizeof(float) * self->half_fft_size + 1);
  memcpy(self->noise_profile, noise_profile,
         sizeof(float) * self->half_fft_size + 1);

  if (self->denoise_parameters->transient_threshold > 1.f) {
    self->transient_detected =
        transient_detector_run(self->transient_detection,
                               self->denoise_parameters->transient_threshold);
  }

  if (self->denoise_parameters->masking_ceiling_limit > 1.f) {
    compute_alpha_and_beta(
        self, self->denoise_parameters->masking_ceiling_limit, 0.f);
  } else {
    memset(self->alpha, 0, self->half_fft_size + 1);
  }

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->noise_profile[k] = self->noise_profile[k] *
                             self->denoise_parameters->noise_rescale *
                             self->alpha[k];
  }

  spectral_smoothing_run(self->spectrum_smoothing,
                         self->denoise_parameters->release_time);

  if (self->transient_detected &&
      self->denoise_parameters->transient_threshold > 1.f) {
    wiener_subtraction(self);
  } else {
    spectral_gating(self);
  }

  memcpy(gain_spectrum, self->gain_spectrum,
         sizeof(float) * self->half_fft_size + 1);
}

static float max_spectral_value(float *spectrum, uint32_t half_fft_size) {
  float max = spectrum[0];
  for (uint32_t k = 1; k <= half_fft_size; k++) {
    max = fmaxf(spectrum[k], max);
  }
  return max;
}

static float min_spectral_value(float *spectrum, uint32_t half_fft_size) {
  float min = spectrum[0];
  for (uint32_t k = 1; k <= half_fft_size; k++) {
    min = fminf(spectrum[k], min);
  }
  return min;
}

static void wiener_subtraction(GainEstimator *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->noise_profile[k] > FLT_MIN) {
      if (self->signal_spectrum[k] > self->noise_profile[k]) {
        self->gain_spectrum[k] =
            (self->signal_spectrum[k] - self->noise_profile[k]) /
            self->signal_spectrum[k];
      } else {
        self->gain_spectrum[k] = 0.f;
      }
    } else {
      self->gain_spectrum[k] = 1.f;
    }
  }
}

static void spectral_gating(GainEstimator *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->noise_profile[k] > FLT_MIN) {
      if (self->signal_spectrum[k] >= self->noise_profile[k]) {
        self->gain_spectrum[k] = 1.f;
      } else {
        self->gain_spectrum[k] = 0.f;
      }
    } else {
      self->gain_spectrum[k] = 1.f;
    }
  }
}

static void compute_alpha_and_beta(GainEstimator *self,
                                   const float masking_ceiling_limit,
                                   const float masking_floor_limit) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->clean_signal_estimation[k] =
        fmaxf(self->signal_spectrum[k] - self->noise_profile[k], FLT_MIN);
  }

  compute_masking_thresholds(self->masking_estimation, self->signal_spectrum,
                             self->masking_thresholds);

  float max_masked_tmp =
      max_spectral_value(self->masking_thresholds, self->half_fft_size);
  float min_masked_tmp =
      min_spectral_value(self->masking_thresholds, self->half_fft_size);

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->masking_thresholds[k] == max_masked_tmp) {
      self->alpha[k] = ALPHA_MIN;
      self->beta[k] = BETA_MIN;
    }
    if (self->masking_thresholds[k] == min_masked_tmp) {
      self->alpha[k] = masking_ceiling_limit;
      self->beta[k] = masking_floor_limit;
    }
    if (self->masking_thresholds[k] < max_masked_tmp &&
        self->masking_thresholds[k] > min_masked_tmp) {
      const float normalized_value =
          (self->masking_thresholds[k] - min_masked_tmp) /
          (max_masked_tmp - min_masked_tmp);

      self->alpha[k] = (1.f - normalized_value) * ALPHA_MIN +
                       normalized_value * masking_ceiling_limit;
      self->beta[k] = (1.f - normalized_value) * BETA_MIN +
                      normalized_value * masking_floor_limit;
    }
  }
}