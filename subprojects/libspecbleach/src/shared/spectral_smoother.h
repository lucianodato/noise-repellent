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

#ifndef SPECTRAL_SMOOTHER_H
#define SPECTRAL_SMOOTHER_H

#include <stdbool.h>
#include <stdint.h>

typedef enum TimeSmoothingType {
  FIXED = 1,
  TRANSIENT_AWARE = 2,
} TimeSmoothingType;

typedef struct TimeSmoothingParameters {
  float smoothing;
  bool transient_protection_enabled;
} TimeSmoothingParameters;

typedef struct SpectralSmoother SpectralSmoother;

SpectralSmoother *spectral_smoothing_initialize(uint32_t fft_size,
                                                TimeSmoothingType type);
void spectral_smoothing_free(SpectralSmoother *self);
bool spectral_smoothing_run(SpectralSmoother *self,
                            TimeSmoothingParameters parameters,
                            float *signal_spectrum);

#endif