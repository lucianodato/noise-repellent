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

#include "louizou_estimator.h"
#include "../../shared/configurations.h"
#include "../../shared/spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N_SMOOTH 0.7f
#define BETA_AT 0.8f
#define GAMMA 0.998f
#define ALPHA_P 0.2f
#define ALPHA_D 0.95f

#define CROSSOVER_POINT1 1000.f
#define CROSSOVER_POINT2 3000.f
#define BAND_1_GAIN 2.0f
#define BAND_2_GAIN 2.0f
#define BAND_3_GAIN 7.0f

typedef struct {
  float *smoothed_power_spectrum;
  float *local_minimum_spectrum;
  float *speech_present_probability_spectrum;
} FrameSpectrum;

static FrameSpectrum *frame_spectrum_initialize(uint32_t frame_size);
static void frame_spectrum_free(FrameSpectrum *self);
static void compute_auto_thresholds(LouizouEstimator *self,
                                    uint32_t sample_rate,
                                    uint32_t noise_spectrum_size,
                                    uint32_t fft_size);

struct LouizouEstimator {
  uint32_t noise_spectrum_size;

  FrameSpectrum *current;
  FrameSpectrum *previous;

  float *auto_thresholds;
  float *previous_noise_spectrum;
  float *time_frequency_smoothing_constant;
  uint32_t *speech_presence_detection;
};

LouizouEstimator *
louizou_estimator_initialize(const uint32_t noise_spectrum_size,
                             const uint32_t sample_rate,
                             const uint32_t fft_size) {
  LouizouEstimator *self =
      (LouizouEstimator *)calloc(1U, sizeof(LouizouEstimator));

  self->noise_spectrum_size = noise_spectrum_size;

  self->auto_thresholds =
      (float *)calloc(self->noise_spectrum_size, sizeof(float));
  self->time_frequency_smoothing_constant =
      (float *)calloc(self->noise_spectrum_size, sizeof(float));
  self->speech_presence_detection =
      (uint32_t *)calloc(self->noise_spectrum_size, sizeof(uint32_t));
  self->previous_noise_spectrum =
      (float *)calloc(self->noise_spectrum_size, sizeof(float));

  compute_auto_thresholds(self, sample_rate, noise_spectrum_size, fft_size);
  self->current = frame_spectrum_initialize(noise_spectrum_size);
  self->previous = frame_spectrum_initialize(noise_spectrum_size);

  return self;
}

void louizou_estimator_free(LouizouEstimator *self) {
  free(self->auto_thresholds);
  free(self->time_frequency_smoothing_constant);
  free(self->speech_presence_detection);
  free(self->previous_noise_spectrum);

  frame_spectrum_free(self->current);
  frame_spectrum_free(self->previous);

  free(self);
}

bool louizou_estimator_run(LouizouEstimator *self, const float *spectrum,
                           float *noise_spectrum) {
  if (!self || !spectrum || !noise_spectrum) {
    return false;
  }

  for (uint32_t k = 0U; k < self->noise_spectrum_size; k++) {
    self->current->smoothed_power_spectrum[k] =
        N_SMOOTH * self->previous->smoothed_power_spectrum[k] +
        (1.f - N_SMOOTH) * spectrum[k];

    if (self->previous->local_minimum_spectrum[k] <
        self->current->smoothed_power_spectrum[k]) {
      self->current->local_minimum_spectrum[k] =
          GAMMA * self->previous->local_minimum_spectrum[k] +
          ((1.f - GAMMA) / (1.f - BETA_AT)) *
              (self->current->smoothed_power_spectrum[k] -
               BETA_AT * self->previous->smoothed_power_spectrum[k]);
    } else {
      self->current->local_minimum_spectrum[k] =
          self->current->smoothed_power_spectrum[k];
    }

    float noisy_speech_ratio = self->current->smoothed_power_spectrum[k] /
                               self->current->local_minimum_spectrum[k];

    if (noisy_speech_ratio > self->auto_thresholds[k]) {
      self->speech_presence_detection[k] = 1U;
    } else {
      self->speech_presence_detection[k] = 0U;
    }

    self->current->speech_present_probability_spectrum[k] =
        ALPHA_P * self->previous->speech_present_probability_spectrum[k] +
        (1.f - ALPHA_P) * (float)self->speech_presence_detection[k];

    self->time_frequency_smoothing_constant[k] =
        ALPHA_D +
        (1.f - ALPHA_D) * self->current->speech_present_probability_spectrum[k];

    noise_spectrum[k] =
        self->time_frequency_smoothing_constant[k] *
            self->previous_noise_spectrum[k] +
        (1.f - self->time_frequency_smoothing_constant[k]) * spectrum[k];
  }

  memcpy(self->previous_noise_spectrum, noise_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->local_minimum_spectrum,
         self->current->local_minimum_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->smoothed_power_spectrum,
         self->current->smoothed_power_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->speech_present_probability_spectrum,
         self->current->speech_present_probability_spectrum,
         sizeof(float) * self->noise_spectrum_size);

  return true;
}

static FrameSpectrum *frame_spectrum_initialize(const uint32_t frame_size) {
  FrameSpectrum *self = (FrameSpectrum *)calloc(1U, sizeof(FrameSpectrum));

  self->smoothed_power_spectrum = (float *)calloc(frame_size, sizeof(float));
  self->local_minimum_spectrum = (float *)calloc(frame_size, sizeof(float));
  self->speech_present_probability_spectrum =
      (float *)calloc(frame_size, sizeof(float));

  return self;
}

static void frame_spectrum_free(FrameSpectrum *self) {
  free(self->smoothed_power_spectrum);
  free(self->local_minimum_spectrum);
  free(self->speech_present_probability_spectrum);

  free(self);
}

static void compute_auto_thresholds(LouizouEstimator *self,
                                    const uint32_t sample_rate,
                                    const uint32_t noise_spectrum_size,
                                    const uint32_t fft_size) {
  int LF = freq_to_fft_bin(CROSSOVER_POINT1, sample_rate, fft_size);
  int MF = freq_to_fft_bin(CROSSOVER_POINT2, sample_rate, fft_size);
  for (int k = 0U; k < noise_spectrum_size; k++) {
    if (k <= LF) {
      self->auto_thresholds[k] = BAND_1_GAIN;
    }
    if (k > LF && k < MF) {
      self->auto_thresholds[k] = BAND_2_GAIN;
    }
    if (k >= MF) {
      self->auto_thresholds[k] = BAND_3_GAIN;
    }
  }
}