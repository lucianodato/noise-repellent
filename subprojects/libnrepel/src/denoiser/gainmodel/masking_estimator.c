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

#include "masking_estimator.h"
#include "../../shared/configurations.h"
#include "../../shared/critical_bands.h"
#include "../../shared/spectral_utils.h"
#include "../../shared/spl_spectrum_converter.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void compute_spectral_spreading_function(MaskingEstimator *self);
static void compute_absolute_thresholds(MaskingEstimator *self);
static float compute_tonality_factor(MaskingEstimator *self,
                                     const float *spectrum, uint32_t band);

struct MaskingEstimator {

  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;

  SplSpectrumConverter *reference_spectrum;
  CriticalBands *critical_bands;

  float *absolute_thresholds;
  float *spectral_spreading_function;
  float *unity_gain_bark_spectrum;
  float *spreaded_unity_gain_bark_spectrum;
  float *threshold_j;
  float *masking_offset;
  float *spreaded_spectrum;
};

MaskingEstimator *masking_estimation_initialize(const uint32_t fft_size,
                                                const uint32_t sample_rate) {

  MaskingEstimator *self =
      (MaskingEstimator *)calloc(1U, sizeof(MaskingEstimator));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2U;
  self->sample_rate = sample_rate;

  self->absolute_thresholds =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));

  self->spectral_spreading_function =
      (float *)calloc((size_t)(N_BARK_BANDS * N_BARK_BANDS), sizeof(float));
  self->unity_gain_bark_spectrum = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->spreaded_unity_gain_bark_spectrum =
      (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->threshold_j = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->masking_offset = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->spreaded_spectrum = (float *)calloc(N_BARK_BANDS, sizeof(float));

  self->reference_spectrum = reference_spectrum_initialize(self->sample_rate);
  self->critical_bands = critical_bands_initialize(
      sample_rate, self->half_fft_size + 1U, N_BARK_BANDS, BARK_SCALE);

  compute_absolute_thresholds(self);
  compute_spectral_spreading_function(self);
  initialize_spectrum_with_value(self->unity_gain_bark_spectrum, N_BARK_BANDS,
                                 1.F);
  direct_matrix_to_vector_spectral_convolution(
      self->spectral_spreading_function, self->unity_gain_bark_spectrum,
      self->spreaded_unity_gain_bark_spectrum, N_BARK_BANDS);

  return self;
}

void masking_estimation_free(MaskingEstimator *self) {
  reference_spectrum_free(self->reference_spectrum);
  critical_bands_free(self->critical_bands);

  free(self->absolute_thresholds);
  free(self->spectral_spreading_function);
  free(self->unity_gain_bark_spectrum);
  free(self->spreaded_unity_gain_bark_spectrum);
  free(self->threshold_j);
  free(self->masking_offset);
  free(self->spreaded_spectrum);

  free(self);
}

bool compute_masking_thresholds(MaskingEstimator *self, const float *spectrum,
                                float *masking_thresholds) {
  if (!self || !spectrum || !masking_thresholds) {
    return false;
  }

  compute_critical_bands_spectrum(self->critical_bands, spectrum);

  direct_matrix_to_vector_spectral_convolution(
      self->spectral_spreading_function,
      get_critical_bands_spectrum(self->critical_bands),
      self->spreaded_spectrum, N_BARK_BANDS);

  for (uint32_t j = 0U; j < N_BARK_BANDS; j++) {

    const float tonality_factor = compute_tonality_factor(self, spectrum, j);

    self->masking_offset[j] = (tonality_factor * (14.5F + (float)(j + 1)) +
                               5.5F * (1.F - tonality_factor));

#if BIAS
    masking_offset[j] = relative_thresholds[j];

    if (j > 15) {
      masking_offset[j] += HIGH_FREQ_BIAS;
    }
#endif

    self->threshold_j[j] =
        powf(10.F, log10f(self->spreaded_spectrum[j]) -
                       (self->masking_offset[j] / 10.F)) -
        (10.F * log10f(self->spreaded_unity_gain_bark_spectrum[j]));

    uint32_t start_pos = 0U;
    if (j == 0) {
      start_pos = 0U;
    } else {
      start_pos = get_intermediate_band_bins(self->critical_bands)[j - 1];
    }

    const uint32_t end_pos =
        get_intermediate_band_bins(self->critical_bands)[j];

    for (uint32_t k = start_pos; k < end_pos; k++) {
      masking_thresholds[k] = self->threshold_j[j];
    }
  }

  convert_spectrum_to_dbspl(self->reference_spectrum, masking_thresholds);

  for (uint32_t k = 1U; k <= self->half_fft_size; k++) {
    masking_thresholds[k] =
        fmaxf(masking_thresholds[k], self->absolute_thresholds[k]);
  }

  return true;
}

static void compute_absolute_thresholds(MaskingEstimator *self) {
  for (uint32_t k = 1U; k <= self->half_fft_size; k++) {

    const float frequency =
        fft_bin_to_freq(k, self->sample_rate, self->half_fft_size);
    self->absolute_thresholds[k] =
        3.64F * powf((frequency / 1000.F), -0.8F) -
        6.5F * expf(-0.6F * powf((frequency / 1000.F - 3.3F), 2.F)) +
        powf(10.F, -3.F) * powf((frequency / 1000.F), 4.F);
  }
}

static void compute_spectral_spreading_function(MaskingEstimator *self) {
  for (uint32_t i = 0U; i < N_BARK_BANDS; i++) {
    for (uint32_t j = 0U; j < N_BARK_BANDS; j++) {
      const uint32_t y = (i + 1) - (j + 1);

      self->spectral_spreading_function[i * N_BARK_BANDS + j] =
          15.81F + 7.5F * ((float)y + 0.474F) -
          17.5F * sqrtf(1.F + ((float)y + 0.474F) * ((float)y + 0.474F));

      self->spectral_spreading_function[i * N_BARK_BANDS + j] = powf(
          10.F, self->spectral_spreading_function[i * N_BARK_BANDS + j] / 10.F);
    }
  }
}

static float compute_tonality_factor(MaskingEstimator *self,
                                     const float *spectrum, uint32_t band) {
  float sum_bins = 0.F;
  float sum_log_bins = 0.F;
  uint32_t start_pos = 0U;
  uint32_t end_pos = 0U;

  if (band == 0) {
    start_pos = band;
    end_pos = get_n_bins_per_band(self->critical_bands)[band];
  } else {
    start_pos = get_intermediate_band_bins(self->critical_bands)[band - 1];
    end_pos = get_intermediate_band_bins(self->critical_bands)[band - 1] +
              get_n_bins_per_band(self->critical_bands)[band];
  }

  for (uint32_t k = start_pos; k < end_pos; k++) {
    sum_bins += spectrum[k];
    sum_log_bins += log10f(spectrum[k]);
  }

  const float SFM =
      10.F *
      (sum_log_bins / (float)(get_n_bins_per_band(self->critical_bands)[band]) -
       log10f(sum_bins /
              (float)(get_n_bins_per_band(self->critical_bands)[band])));

  const float tonality_factor = fminf(SFM / -60.F, 1.F);

  return tonality_factor;
}
