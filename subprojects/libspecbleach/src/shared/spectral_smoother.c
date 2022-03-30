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

static void spectrum_time_smoothing(SpectralSmoother *self, float smoothing);
static void spectrum_transient_aware_time_smoothing(SpectralSmoother *self,
                                                    float smoothing,
                                                    float *spectrum);

struct SpectralSmoother {
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  float adaptive_coefficient;
  float previous_adaptive_coefficient;

  float *noise_spectrum;

  float *smoothed_spectrum;
  float *smoothed_spectrum_previous;

  TimeSmoothingType type;
  TransientDetector *transient_detection;
};

SpectralSmoother *spectral_smoothing_initialize(const uint32_t fft_size,
                                                TimeSmoothingType type) {
  SpectralSmoother *self =
      (SpectralSmoother *)calloc(1U, sizeof(SpectralSmoother));

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->type = type;
  self->previous_adaptive_coefficient = 0.F;
  self->adaptive_coefficient = 0.F;

  self->noise_spectrum =
      (float *)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_spectrum =
      (float *)calloc(self->real_spectrum_size, sizeof(float));
  self->smoothed_spectrum_previous =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

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

bool spectral_smoothing_run(SpectralSmoother *self,
                            TimeSmoothingParameters parameters,
                            float *signal_spectrum) {
  if (!self || !signal_spectrum) {
    return false;
  }

  memcpy(self->smoothed_spectrum, signal_spectrum,
         sizeof(float) * self->real_spectrum_size);

  switch (self->type) {
  case FIXED_RELEASE:
    spectrum_time_smoothing(self, parameters.smoothing);
    break;
  case TRANSIENT_AWARE:
    if (parameters.transient_protection_enabled) {
      spectrum_transient_aware_time_smoothing(self, parameters.smoothing,
                                              signal_spectrum);
    } else {
      spectrum_time_smoothing(self, parameters.smoothing);
    }
    break;
  default:
    break;
  }

  memcpy(self->smoothed_spectrum_previous, self->smoothed_spectrum,
         sizeof(float) * self->real_spectrum_size);
  memcpy(signal_spectrum, self->smoothed_spectrum,
         sizeof(float) * self->real_spectrum_size);

  return true;
}

static void spectrum_transient_aware_time_smoothing(SpectralSmoother *self,
                                                    const float smoothing,
                                                    float *spectrum) {

  if (!transient_detector_run(self->transient_detection, spectrum)) {
    for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
      if (self->smoothed_spectrum[k] > self->smoothed_spectrum_previous[k]) {
        self->smoothed_spectrum[k] =
            smoothing * self->smoothed_spectrum_previous[k] +
            (1.F - smoothing) * self->smoothed_spectrum[k];
      }
    }
  }
}

static void spectrum_time_smoothing(SpectralSmoother *self,
                                    const float smoothing) {
  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
    if (self->smoothed_spectrum[k] > self->smoothed_spectrum_previous[k]) {
      self->smoothed_spectrum[k] =
          smoothing * self->smoothed_spectrum_previous[k] +
          (1.F - smoothing) * self->smoothed_spectrum[k];
    }
  }
}