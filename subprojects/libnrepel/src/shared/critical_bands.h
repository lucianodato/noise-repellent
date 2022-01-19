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

#ifndef CRITICAL_BANDS_H
#define CRITICAL_BANDS_H

#include <stdbool.h>
#include <stdint.h>

typedef struct CriticalBands CriticalBands;

typedef enum CriticalBandType {
  BARK_SCALE = 0,
  MEL_SCALE = 1,
  ERB_SCALE = 2,
  // TODO (luciano/todo): opus critical bands
  // TODO (luciano/todo): octave bands
} CriticalBandType;

typedef struct CriticalBandIndexes {
  uint32_t start_position;
  uint32_t end_position;
} CriticalBandIndexes;

CriticalBands *critical_bands_initialize(uint32_t sample_rate,
                                         uint32_t spectrum_size,
                                         uint32_t number_bands,
                                         CriticalBandType type);
void critical_bands_free(CriticalBands *self);
bool compute_critical_bands_spectrum(CriticalBands *self, const float *spectrum,
                                     float *critical_bands);
CriticalBandIndexes get_band_indexes(CriticalBands *self, uint32_t band_number);

#endif