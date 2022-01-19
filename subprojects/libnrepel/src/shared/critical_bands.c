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
#include <math.h>
#include <stdlib.h>

static void compute_mapping_spectrum(CriticalBands *self);
static void compute_band_indexes(CriticalBands *self);

struct CriticalBands {
  float *mapping_spectrum;
  uint32_t *intermediate_band_bins;
  uint32_t *n_bins_per_band;

  uint32_t spectrum_size;
  uint32_t sample_rate;
  uint32_t number_bands;
  CriticalBandType type;
  CriticalBandIndexes band_indexes;
};

CriticalBands *critical_bands_initialize(const uint32_t sample_rate,
                                         const uint32_t spectrum_size,
                                         const uint32_t number_bands,
                                         const CriticalBandType type) {

  CriticalBands *self =
      (CriticalBands *)calloc(number_bands, sizeof(CriticalBands));

  self->spectrum_size = spectrum_size;
  self->number_bands = number_bands;
  self->sample_rate = sample_rate;
  self->type = type;

  self->mapping_spectrum =
      (float *)calloc(self->spectrum_size + 1U, sizeof(float));
  self->intermediate_band_bins =
      (uint32_t *)calloc(self->number_bands, sizeof(uint32_t));
  self->n_bins_per_band =
      (uint32_t *)calloc(self->number_bands, sizeof(uint32_t));

  compute_mapping_spectrum(self);
  compute_band_indexes(self);

  return self;
};

void critical_bands_free(CriticalBands *self) {
  free(self->intermediate_band_bins);
  free(self->n_bins_per_band);
  free(self->mapping_spectrum);

  free(self);
}

static void compute_mapping_spectrum(CriticalBands *self) {
  for (uint32_t k = 0U; k <= self->spectrum_size; k++) {
    const float frequency =
        fft_bin_to_freq(k, self->sample_rate, self->spectrum_size);
    switch (self->type) {
    case BARK_SCALE:
      self->mapping_spectrum[k] =
          floorf(13.F * atanf(0.00076F * frequency) +
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

    uint32_t current_band = j + 1U;
    while (self->mapping_spectrum[last_position + counter] <=
               (float)current_band &&
           (last_position + counter) <= self->spectrum_size) {
      counter++;
    }

    self->n_bins_per_band[j] = counter;
    self->intermediate_band_bins[j] = last_position;

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
  uint32_t start_pos = 0U;
  uint32_t end_pos = 0U;

  if (band_number == 0) {
    start_pos = 1U;
    end_pos = self->n_bins_per_band[band_number];
  } else {
    start_pos = self->intermediate_band_bins[band_number - 1];
    end_pos = self->intermediate_band_bins[band_number - 1] +
              self->n_bins_per_band[band_number];
  }

  return (CriticalBandIndexes){
      .start_position = start_pos,
      .end_position = end_pos,
  };
}