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
#include "../shared/configurations.h"
#include "../shared/general_utils.h"
#include "../shared/spectral_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct FrameSpectrum {
  float *smoothed_spectrum;
  float *local_minimum_spectrum;
  float *speech_present_probability_spectrum;
} FrameSpectrum;

static FrameSpectrum *frame_spectrum_initialize(uint32_t frame_size);
static void frame_spectrum_free(FrameSpectrum *self);
static void compute_auto_thresholds(LouizouEstimator *self,
                                    uint32_t sample_rate,
                                    uint32_t noise_spectrum_size,
                                    uint32_t fft_size);
static void update_frame_spectums(LouizouEstimator *self,
                                  const float *noise_spectrum);

struct LouizouEstimator {
  uint32_t noise_spectrum_size;
  float noisy_speech_ratio;

  FrameSpectrum *current;
  FrameSpectrum *previous;

  float *minimum_detection_thresholds;
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

  self->minimum_detection_thresholds =
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

  self->noisy_speech_ratio = 0.F;

  return self;
}

void louizou_estimator_free(LouizouEstimator *self) {
  free(self->minimum_detection_thresholds);
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

  for (uint32_t k = 1U; k < self->noise_spectrum_size; k++) {
    self->current->smoothed_spectrum[k] =
        N_SMOOTH * self->previous->smoothed_spectrum[k] +
        (1.F - N_SMOOTH) * spectrum[k];

    if (self->previous->local_minimum_spectrum[k] <
        self->current->smoothed_spectrum[k]) {
      self->current->local_minimum_spectrum[k] =
          GAMMA * self->previous->local_minimum_spectrum[k] +
          ((1.F - GAMMA) / (1.F - BETA_AT)) *
              (self->current->smoothed_spectrum[k] -
               BETA_AT * self->previous->smoothed_spectrum[k]);
    } else {
      self->current->local_minimum_spectrum[k] =
          self->current->smoothed_spectrum[k];
    }

    self->noisy_speech_ratio =
        sanitize_denormal(self->current->smoothed_spectrum[k] /
                          self->current->local_minimum_spectrum[k]);

    if (self->noisy_speech_ratio > self->minimum_detection_thresholds[k]) {
      self->speech_presence_detection[k] = 1U;
    } else {
      self->speech_presence_detection[k] = 0U;
    }

    self->current->speech_present_probability_spectrum[k] =
        ALPHA_P * self->previous->speech_present_probability_spectrum[k] +
        (1.F - ALPHA_P) * (float)self->speech_presence_detection[k];

    self->time_frequency_smoothing_constant[k] =
        ALPHA_D +
        (1.F - ALPHA_D) * self->current->speech_present_probability_spectrum[k];

    noise_spectrum[k] =
        self->time_frequency_smoothing_constant[k] *
            self->previous_noise_spectrum[k] +
        (1.F - self->time_frequency_smoothing_constant[k]) * spectrum[k];
  }

  update_frame_spectums(self, noise_spectrum);

  return true;
}

static void update_frame_spectums(LouizouEstimator *self,
                                  const float *noise_spectrum) {
  memcpy(self->previous_noise_spectrum, noise_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->local_minimum_spectrum,
         self->current->local_minimum_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->smoothed_spectrum, self->current->smoothed_spectrum,
         sizeof(float) * self->noise_spectrum_size);
  memcpy(self->previous->speech_present_probability_spectrum,
         self->current->speech_present_probability_spectrum,
         sizeof(float) * self->noise_spectrum_size);
}

static FrameSpectrum *frame_spectrum_initialize(const uint32_t frame_size) {
  FrameSpectrum *self = (FrameSpectrum *)calloc(1U, sizeof(FrameSpectrum));

  self->smoothed_spectrum = (float *)calloc(frame_size, sizeof(float));
  self->local_minimum_spectrum = (float *)calloc(frame_size, sizeof(float));
  self->speech_present_probability_spectrum =
      (float *)calloc(frame_size, sizeof(float));

  initialize_spectrum_with_value(self->local_minimum_spectrum, frame_size,
                                 FLT_MIN);

  return self;
}

static void frame_spectrum_free(FrameSpectrum *self) {
  free(self->smoothed_spectrum);
  free(self->local_minimum_spectrum);
  free(self->speech_present_probability_spectrum);

  free(self);
}

static void compute_auto_thresholds(LouizouEstimator *self,
                                    const uint32_t sample_rate,
                                    const uint32_t noise_spectrum_size,
                                    const uint32_t fft_size) {
  uint32_t LF = freq_to_fft_bin(CROSSOVER_POINT1, sample_rate, fft_size);
  uint32_t MF = freq_to_fft_bin(CROSSOVER_POINT2, sample_rate, fft_size);
  for (uint32_t k = 0U; k < noise_spectrum_size; k++) {
    if (k <= LF) {
      self->minimum_detection_thresholds[k] = BAND_1_LEVEL;
    }
    if (k > LF && k < MF) {
      self->minimum_detection_thresholds[k] = BAND_2_LEVEL;
    }
    if (k >= MF) {
      self->minimum_detection_thresholds[k] = BAND_3_LEVEL;
    }
  }
}