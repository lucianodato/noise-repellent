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

#ifndef MASKING_ESTIMATOR_H
#define MASKING_ESTIMATOR_H

#include "spectral_features.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct MaskingEstimator MaskingEstimator;

MaskingEstimator *masking_estimation_initialize(uint32_t fft_size,
                                                uint32_t sample_rate,
                                                SpectrumType spectrum_type);
void masking_estimation_free(MaskingEstimator *self);
bool compute_masking_thresholds(MaskingEstimator *self, const float *spectrum,
                                float *masking_thresholds);

#endif