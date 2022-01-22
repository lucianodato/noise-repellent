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

#ifndef ABSOLUTE_HEARING_THRESHOLDS_H
#define ABSOLUTE_HEARING_THRESHOLDS_H

#include "spectral_features.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct AbsoluteHearingThresholds AbsoluteHearingThresholds;

AbsoluteHearingThresholds *
absolute_hearing_thresholds_initialize(uint32_t sample_rate, uint32_t fft_size,
                                       SpectrumType spectrum_type);
void absolute_hearing_thresholds_free(AbsoluteHearingThresholds *self);
bool apply_thresholds_as_floor(AbsoluteHearingThresholds *self,
                               float *spectrum);

#endif