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

#include "absolute_hearing_thresholds.h"
#include "configurations.h"
#include "fft_transform.h"
#include "spectral_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void generate_sinewave(AbsoluteHearingThresholds *self);
static void compute_spl_reference_spectrum(AbsoluteHearingThresholds *self);
static void compute_absolute_thresholds(AbsoluteHearingThresholds *self);

struct AbsoluteHearingThresholds {
  float *sinewave;
  float *window;
  float *spl_reference_values;
  float *absolute_thresholds;

  SpectralFeatures *spectral_features;
  FftTransform *fft_transform;
  SpectrumType spectrum_type;

  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  float sine_wave_amplitude;
  float sine_wave_frequency;
  float reference_level;
};

AbsoluteHearingThresholds *
absolute_hearing_thresholds_initialize(const uint32_t sample_rate,
                                       const uint32_t fft_size,
                                       SpectrumType spectrum_type) {
  AbsoluteHearingThresholds *self = (AbsoluteHearingThresholds *)calloc(
      1U, sizeof(AbsoluteHearingThresholds));

  self->fft_transform = fft_transform_initialize_bins(fft_size);

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->sample_rate = sample_rate;
  self->spectrum_type = spectrum_type;
  self->sine_wave_amplitude = SINE_AMPLITUDE;
  self->sine_wave_frequency = REFERENCE_SINE_WAVE_FREQ;
  self->reference_level = REFERENCE_LEVEL;

  self->spl_reference_values =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

  self->absolute_thresholds =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

  self->sinewave = (float *)calloc(self->fft_size, sizeof(float));
  self->window = (float *)calloc(self->fft_size, sizeof(float));

  self->spectral_features =
      spectral_features_initialize(self->real_spectrum_size);

  generate_sinewave(self);
  get_fft_window(self->window, self->fft_size, VORBIS_WINDOW);
  compute_spl_reference_spectrum(self);
  compute_absolute_thresholds(self);

  return self;
}

void absolute_hearing_thresholds_free(AbsoluteHearingThresholds *self) {
  fft_transform_free(self->fft_transform);
  spectral_features_free(self->spectral_features);

  free(self->sinewave);
  free(self->window);
  free(self->spl_reference_values);
  free(self->absolute_thresholds);

  free(self);
}

static void generate_sinewave(AbsoluteHearingThresholds *self) {
  for (uint32_t k = 0U; k < self->fft_size; k++) {
    self->sinewave[k] =
        self->sine_wave_amplitude *
        sinf((2.F * M_PI * (float)k * self->sine_wave_frequency) /
             (float)self->sample_rate);
  }
}

static void compute_spl_reference_spectrum(AbsoluteHearingThresholds *self) {
  for (uint32_t k = 0U; k < self->fft_size; k++) {
    get_fft_input_buffer(self->fft_transform)[k] =
        self->sinewave[k] * self->window[k];
  }

  compute_forward_fft(self->fft_transform);

  float *reference_spectrum = get_spectral_feature(
      self->spectral_features, get_fft_output_buffer(self->fft_transform),
      self->fft_size, self->spectrum_type);

  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
    self->spl_reference_values[k] =
        self->reference_level - 10.F * log10f(reference_spectrum[k]);
  }
}

bool apply_thresholds_as_floor(AbsoluteHearingThresholds *self,
                               float *spectrum) {
  if (!self || !spectrum) {
    return false;
  }

  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
    spectrum[k] = fmaxf(spectrum[k] + self->spl_reference_values[k],
                        self->absolute_thresholds[k]);
  }

  return true;
}

static void compute_absolute_thresholds(AbsoluteHearingThresholds *self) {
  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {

    const float frequency =
        fft_bin_to_freq(k, self->sample_rate, self->fft_size);
    self->absolute_thresholds[k] =
        3.64F * powf((frequency / 1000.F), -0.8F) -
        6.5F * expf(-0.6F * powf((frequency / 1000.F - 3.3F), 2.F)) +
        powf(10.F, -3.F) * powf((frequency / 1000.F), 4.F);
  }
}
