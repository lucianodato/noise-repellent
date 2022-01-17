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

#ifndef FFT_TRANSFORM_H
#define FFT_TRANSFORM_H

#include <stdbool.h>
#include <stdint.h>

typedef struct CriticalBands CriticalBands;

typedef enum CriticalBandType {
  BARK_SCALE = 0,
  MEL_SCALE = 1,
  ERB_SCALE = 2,
} CriticalBandType;

CriticalBands *critical_bands_initialize(uint32_t sample_rate,
                                         uint32_t frame_size,
                                         uint32_t number_bands,
                                         CriticalBandType type);
void critical_bands_free(CriticalBands *self);
bool compute_critical_bands_spectrum(CriticalBands *self,
                                     const float *spectrum);
float *get_critical_bands_spectrum(CriticalBands *self);
uint32_t *get_intermediate_band_bins(CriticalBands *self);
uint32_t *get_n_bins_per_band(CriticalBands *self);

#endif