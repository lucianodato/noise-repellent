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

#include "spectral_smoother.h"
#include "transient_detector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void get_release_coefficient(SpectralSmoother *self, float release);
static void spectrum_time_smoothing(SpectralSmoother *self);
static void spectrum_transient_aware_time_smoothing(SpectralSmoother *self,
                                                    float *spectrum);
static void spectrum_adaptive_time_smoothing(SpectralSmoother *self,
                                             float *spectrum,
                                             const float *noise_thresholds);

struct SpectralSmoother {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  uint32_t hop;
  float adaptive_coefficient;
  float previous_adaptive_coefficient;

  float *noise_spectrum;

  float *smoothed_spectrum;
  float *smoothed_spectrum_previous;

  float release_coefficient;
  TimeSmoothingType type;
  TransientDetector *transient_detection;
};

SpectralSmoother *spectral_smoothing_initialize(const uint32_t fft_size,
                                                const uint32_t sample_rate,
                                                const uint32_t hop,
                                                TimeSmoothingType type) {
  SpectralSmoother *self =
      (SpectralSmoother *)calloc(1U, sizeof(SpectralSmoother));

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->sample_rate = sample_rate;
  self->hop = hop;
  self->type = type;
  self->previous_adaptive_coefficient = 0.F;
  self->adaptive_coefficient = 0.F;

  self->noise_spectrum =
      (float *)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_spectrum =
      (float *)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_spectrum_previous =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

  self->release_coefficient = 0.F;

  self->transient_detection = transient_detector_initialize(self->fft_size);

  return self;
}

void spectral_smoothing_free(SpectralSmoother *self) {
  free(self->noise_spectrum);
  free(self->smoothed_spectrum);
  free(self->smoothed_spectrum_previous);

  transient_detector_free(self->transient_detection);

  free(self);
}

bool spectral_smoothing_run(SpectralSmoother *self, const float release,
                            float *signal_spectrum,
                            const float *noise_spectrum) {
  if (!self || !signal_spectrum) {
    return false;
  }
  get_release_coefficient(self, release);

  memcpy(self->smoothed_spectrum, signal_spectrum,
         sizeof(float) * self->real_spectrum_size);

  switch (self->type) {
  case FIXED_RELEASE:
    spectrum_time_smoothing(self);
    break;
  case ADAPTIVE_RELEASE:
    spectrum_adaptive_time_smoothing(self, signal_spectrum, noise_spectrum);
    break;
  case TRANSIENT_AWARE:
    spectrum_transient_aware_time_smoothing(self, signal_spectrum);
  default:
    break;
  }

  memcpy(self->smoothed_spectrum_previous, self->smoothed_spectrum,
         sizeof(float) * self->real_spectrum_size);
  memcpy(signal_spectrum, self->smoothed_spectrum,
         sizeof(float) * self->real_spectrum_size);

  return true;
}

static void get_release_coefficient(SpectralSmoother *self,
                                    const float release) {
  if (release != 0.F) {
    self->release_coefficient = expf(
        -1000.F / (((release) * (float)self->sample_rate) / (float)self->hop));
  } else {
    self->release_coefficient = 0.F;
  }
}

static void spectrum_transient_aware_time_smoothing(SpectralSmoother *self,
                                                    float *spectrum) {

  if (!transient_detector_run(self->transient_detection, spectrum)) {
    for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
      if (self->smoothed_spectrum[k] > self->smoothed_spectrum_previous[k]) {
        self->smoothed_spectrum[k] =
            self->release_coefficient * self->smoothed_spectrum_previous[k] +
            (1.F - self->release_coefficient) * self->smoothed_spectrum[k];
      }
    }
  }
}

static void spectrum_time_smoothing(SpectralSmoother *self) {
  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
    if (self->smoothed_spectrum[k] > self->smoothed_spectrum_previous[k]) {
      self->smoothed_spectrum[k] =
          self->release_coefficient * self->smoothed_spectrum_previous[k] +
          (1.F - self->release_coefficient) * self->smoothed_spectrum[k];
    }
  }
}

static void spectrum_adaptive_time_smoothing(SpectralSmoother *self,
                                             float *spectrum,
                                             const float *noise_thresholds) {
  float numerator = 0.F;
  float denominator = 0.F;
  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
    numerator += fabsf(spectrum[k] - noise_thresholds[k]);
    denominator += noise_thresholds[k];
  }
  float spectrum_discrepancy = numerator / denominator;
  self->adaptive_coefficient = fminf(spectrum_discrepancy, 1.F);

  float current_coefficient = 0.F;
  if (self->adaptive_coefficient > self->previous_adaptive_coefficient) {
    current_coefficient = self->release_coefficient;
  }

  float smoothed_coefficient =
      current_coefficient * self->previous_adaptive_coefficient +
      (1.F - current_coefficient) * self->adaptive_coefficient;

  if (self->release_coefficient == 0.F) {
    smoothed_coefficient = 0.F;
  }

  self->previous_adaptive_coefficient = smoothed_coefficient;

  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
    self->smoothed_spectrum[k] =
        (1.F - smoothed_coefficient) * self->smoothed_spectrum[k] +
        smoothed_coefficient * self->smoothed_spectrum_previous[k];
  }
}