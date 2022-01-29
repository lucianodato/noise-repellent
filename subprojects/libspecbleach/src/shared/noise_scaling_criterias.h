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

#ifndef OVERSUBTRACTION_CRITERIAS_H
#define OVERSUBTRACTION_CRITERIAS_H

#include "critical_bands.h"
#include "spectral_features.h"
#include <stdbool.h>
#include <stdint.h>

typedef enum NoiseScalingType {
  A_POSTERIORI_SNR_CRITICAL_BANDS = 0,
  A_POSTERIORI_SNR = 1,
  MASKING_THRESHOLDS = 2,
} NoiseScalingType;

typedef struct NoiseScalingParameters {
  float undersubtraction;
  float oversubtraction;
} NoiseScalingParameters;

typedef struct NoiseScalingCriterias NoiseScalingCriterias;

NoiseScalingCriterias *noise_scaling_criterias_initialize(
    NoiseScalingType subtraction_type, uint32_t fft_size,
    CriticalBandType critical_band_type, uint32_t sample_rate,
    SpectrumType spectrum_type);
void noise_scaling_criterias_free(NoiseScalingCriterias *self);
bool apply_noise_scaling_criteria(NoiseScalingCriterias *self,
                                  const float *spectrum, float *noise_spectrum,
                                  float *alpha, float *beta,
                                  NoiseScalingParameters parameters);

#endif