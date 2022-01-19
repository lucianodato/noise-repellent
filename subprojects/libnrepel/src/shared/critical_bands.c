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

#include "critical_bands.h"
#include "configurations.h"
#include "spectral_utils.h"
#include <math.h>
#include <stdlib.h>

static void compute_mapping_spectrum(CriticalBands *self);
static void compute_band_indexes(CriticalBands *self);

struct CriticalBands {
  float *mapping_spectrum;
  uint32_t *band_delimiter_bins;
  uint32_t *number_bins_per_band;

  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;
  uint32_t number_bands;
  CriticalBandType type;
  CriticalBandIndexes band_indexes;
};

CriticalBands *critical_bands_initialize(const uint32_t sample_rate,
                                         const uint32_t fft_size,
                                         const uint32_t number_bands,
                                         const CriticalBandType type) {

  CriticalBands *self =
      (CriticalBands *)calloc(number_bands, sizeof(CriticalBands));

  self->fft_size = fft_size;
  self->half_fft_size = fft_size / 2U;
  self->number_bands = number_bands;
  self->sample_rate = sample_rate;
  self->type = type;

  self->mapping_spectrum =
      (float *)calloc(self->half_fft_size + 1U, sizeof(float));
  self->band_delimiter_bins =
      (uint32_t *)calloc(self->number_bands, sizeof(uint32_t));
  self->number_bins_per_band =
      (uint32_t *)calloc(self->number_bands, sizeof(uint32_t));

  compute_mapping_spectrum(self);
  compute_band_indexes(self);

  return self;
};

void critical_bands_free(CriticalBands *self) {
  free(self->band_delimiter_bins);
  free(self->number_bins_per_band);
  free(self->mapping_spectrum);

  free(self);
}

static void compute_mapping_spectrum(CriticalBands *self) {
  for (uint32_t k = 1U; k <= self->half_fft_size; k++) {
    const float frequency =
        fft_bin_to_freq(k, self->sample_rate, self->fft_size);
    switch (self->type) {
    case BARK_SCALE:
      self->mapping_spectrum[k] =
          truncf(13.F * atanf(0.00076F * frequency) +
                 3.5F * atanf(powf(frequency / 7500.F, 2.F)));
      break;
    case MEL_SCALE: // FIXME (luciano/todo): broken mapping
      self->mapping_spectrum[k] = 2595.F * log10f(1.F + frequency / 700.F);
      break;
    case ERB_SCALE: // FIXME (luciano/todo): broken mapping
      self->mapping_spectrum[k] =
          6.23F * powf(frequency, 2.F) + 93.39F * frequency + 28.52F;
      break;
    default:
      break;
    }
  }
}

static void compute_band_indexes(CriticalBands *self) {

  uint32_t last_position = 1U;
  for (uint32_t j = 0U; j < self->number_bands; j++) {
    uint32_t counter = 0U;

    float current_band = (float)(j + 1U);
    while (self->mapping_spectrum[last_position + counter] == current_band) {
      counter++;
    }

    self->number_bins_per_band[j] = counter;
    self->band_delimiter_bins[j] = last_position;

    last_position += counter;
  }
}

bool compute_critical_bands_spectrum(CriticalBands *self, const float *spectrum,
                                     float *critical_bands) {
  if (!spectrum) {
    return false;
  }

  for (uint32_t j = 0U; j < self->number_bands; j++) {

    self->band_indexes = get_band_indexes(self, j);

    for (uint32_t k = self->band_indexes.start_position;
         k < self->band_indexes.end_position; k++) {
      critical_bands[j] += spectrum[k];
    }
  }

  return true;
}

CriticalBandIndexes get_band_indexes(CriticalBands *self,
                                     const uint32_t band_number) {
  return (CriticalBandIndexes){
      .start_position = self->band_delimiter_bins[band_number],
      .end_position = self->band_delimiter_bins[band_number] +
                      self->number_bins_per_band[band_number],
  };
}