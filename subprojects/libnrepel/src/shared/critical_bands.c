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

struct CriticalBands {
  float *mapping_spectrum;
  float *converted_spectrum;
  uint32_t *intermediate_band_bins;
  uint32_t *n_bins_per_band;

  uint32_t frame_size;
  uint32_t sample_rate;
  uint32_t number_bands;
  CriticalBandType type;
};

CriticalBands *critical_bands_initialize(const uint32_t sample_rate,
                                         const uint32_t frame_size,
                                         const uint32_t number_bands,
                                         const CriticalBandType type) {

  CriticalBands *self =
      (CriticalBands *)calloc(number_bands, sizeof(CriticalBands));

  self->frame_size = frame_size;
  self->number_bands = number_bands;
  self->sample_rate = sample_rate;
  self->type = type;

  self->mapping_spectrum = (float *)calloc(frame_size + 1U, sizeof(float));
  self->converted_spectrum = (float *)calloc(number_bands, sizeof(float));
  self->intermediate_band_bins =
      (uint32_t *)calloc(number_bands, sizeof(uint32_t));
  self->n_bins_per_band = (uint32_t *)calloc(number_bands, sizeof(uint32_t));

  compute_mapping_spectrum(self);

  return self;
};

void critical_bands_free(CriticalBands *self) {
  free(self->intermediate_band_bins);
  free(self->n_bins_per_band);
  free(self->converted_spectrum);
  free(self->mapping_spectrum);
}

static void compute_mapping_spectrum(CriticalBands *self) {
  for (uint32_t k = 1U; k <= self->frame_size; k++) {
    const float frequency =
        fft_bin_to_freq(k, self->sample_rate, self->frame_size);
    switch (self->type) {
    case BARK_SCALE:
      self->mapping_spectrum[k] = 13.F * atanf(0.00076F * frequency) +
                                  3.5F * atanf(powf(frequency / 7500.F, 2.F));
      break;
    case MEL_SCALE:
      self->mapping_spectrum[k] = 2595.F * log10f(1.F + frequency / 700.F);
      break;
    case ERB_SCALE:
      self->mapping_spectrum[k] =
          6.23F * powf(frequency, 2.F) + 93.39F * frequency + 28.52F;
      break;
    default:
      break;
    }
  }
}

bool compute_critical_bands_spectrum(CriticalBands *self,
                                     const float *spectrum) {
  if (!spectrum) {
    return false;
  }

  uint32_t last_position = 0U;

  for (uint32_t j = 0U; j < N_BARK_BANDS; j++) {
    uint32_t counter = 0U;
    if (j == 0) {
      counter = 1U;
    }

    self->converted_spectrum[j] = 0.F;

    while (floorf(self->mapping_spectrum[last_position + counter]) ==
           ((float)j + 1)) {
      self->converted_spectrum[j] += spectrum[last_position + counter];
      counter++;
    }

    last_position += counter;

    self->n_bins_per_band[j] = counter;
    self->intermediate_band_bins[j] = last_position;
  }

  return true;
}

float *get_critical_bands_spectrum(CriticalBands *self) {
  return self->converted_spectrum;
}
uint32_t *get_intermediate_band_bins(CriticalBands *self) {
  return self->intermediate_band_bins;
}
uint32_t *get_n_bins_per_band(CriticalBands *self) {
  return self->n_bins_per_band;
}