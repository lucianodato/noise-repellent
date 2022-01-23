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

#include "oversubtraction_criterias.h"
#include "configurations.h"
#include "masking_estimator.h"
#include "spectral_utils.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

static void
a_posteriori_snr_critical_bands(OversubtractionCriterias *self,
                                const float *spectrum, float *noise_spectrum,
                                OversustractionParameters parameters);
static void a_posteriori_snr(OversubtractionCriterias *self,
                             const float *spectrum, float *noise_spectrum,
                             OversustractionParameters parameters);
static void masking_thresholds(OversubtractionCriterias *self,
                               const float *spectrum, float *noise_spectrum,
                               OversustractionParameters parameters);

struct OversubtractionCriterias {
  OversubtractionType oversubtraction_type;
  uint32_t fft_size;
  uint32_t real_spectrum_size;
  uint32_t sample_rate;
  SpectrumType spectrum_type;
  uint32_t number_critical_bands;

  float *alpha;
  float *beta;
  float *masking_thresholds;
  float *clean_signal_estimation;

  MaskingEstimator *masking_estimation;
  CriticalBands *critical_bands;
  CriticalBandIndexes band_indexes;
  CriticalBandType critical_band_type;
  float *bark_noise_profile;
  float *bark_reference_spectrum;
};

OversubtractionCriterias *oversubtraction_criterias_initialize(
    const OversubtractionType subtraction_type, const uint32_t fft_size,
    const CriticalBandType critical_band_type, const uint32_t sample_rate,
    SpectrumType spectrum_type) {

  OversubtractionCriterias *self =
      (OversubtractionCriterias *)calloc(1U, sizeof(OversubtractionCriterias));

  self->oversubtraction_type = subtraction_type;
  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->critical_band_type = critical_band_type;
  self->sample_rate = sample_rate;
  self->spectrum_type = spectrum_type;

  self->critical_bands = critical_bands_initialize(
      self->sample_rate, self->fft_size, self->critical_band_type);
  self->masking_estimation = masking_estimation_initialize(
      self->fft_size, self->sample_rate, self->spectrum_type);
  self->number_critical_bands =
      get_number_of_critical_bands(self->critical_bands);

  self->bark_noise_profile =
      (float *)calloc(self->number_critical_bands, sizeof(float));
  self->bark_reference_spectrum =
      (float *)calloc(self->number_critical_bands, sizeof(float));

  self->alpha = (float *)calloc(self->real_spectrum_size, sizeof(float));
  self->beta = (float *)calloc(self->real_spectrum_size, sizeof(float));
  self->masking_thresholds =
      (float *)calloc(self->real_spectrum_size, sizeof(float));
  self->clean_signal_estimation =
      (float *)calloc(self->real_spectrum_size, sizeof(float));

  return self;
}

void oversubtraction_criterias_free(OversubtractionCriterias *self) {
  critical_bands_free(self->critical_bands);
  masking_estimation_free(self->masking_estimation);

  free(self->alpha);
  free(self->beta);
  free(self->clean_signal_estimation);
  free(self->masking_thresholds);
  free(self->bark_noise_profile);
  free(self->bark_reference_spectrum);

  free(self);
}

bool apply_oversustraction_criteria(OversubtractionCriterias *self,
                                    const float *spectrum,
                                    float *noise_spectrum,
                                    OversustractionParameters parameters) {
  if (!spectrum || !noise_spectrum) {
    return false;
  }

  switch (self->oversubtraction_type) {
  case A_POSTERIORI_SNR_CRITICAL_BANDS:
    a_posteriori_snr_critical_bands(self, spectrum, noise_spectrum, parameters);
    break;
  case A_POSTERIORI_SNR:
    a_posteriori_snr(self, spectrum, noise_spectrum, parameters);
    break;
  case MASKING_THRESHOLDS:
    masking_thresholds(self, spectrum, noise_spectrum, parameters);
    break;

  default:
    break;
  }

  return true;
}

static void
a_posteriori_snr_critical_bands(OversubtractionCriterias *self,
                                const float *spectrum, float *noise_spectrum,
                                OversustractionParameters parameters) {

  compute_critical_bands_spectrum(self->critical_bands, noise_spectrum,
                                  self->bark_noise_profile);
  compute_critical_bands_spectrum(self->critical_bands, spectrum,
                                  self->bark_reference_spectrum);

  float a_posteriori_snr = 20.F;
  float oversustraction_factor = 1.F;

  for (uint32_t j = 0U; j < self->number_critical_bands; j++) {

    self->band_indexes = get_band_indexes(self->critical_bands, j);

    a_posteriori_snr = 10.F * log10f(self->bark_reference_spectrum[j] /
                                     self->bark_noise_profile[j]);

    if (a_posteriori_snr >= 0.F && a_posteriori_snr <= 20.F) {
      oversustraction_factor = -0.05F * (a_posteriori_snr) + 5.F;
    } else if (a_posteriori_snr < 0.F) {
      oversustraction_factor = 5.F;
    } else if (a_posteriori_snr > 20.F) {
      oversustraction_factor = 1.F;
    }

    for (uint32_t k = self->band_indexes.start_position;
         k < self->band_indexes.end_position; k++) {
      noise_spectrum[k] *= parameters.noise_rescale * oversustraction_factor;
    }
  }
}

static void a_posteriori_snr(OversubtractionCriterias *self,
                             const float *spectrum, float *noise_spectrum,
                             OversustractionParameters parameters) {
  float a_posteriori_snr = 20.F;
  float oversustraction_factor = 1.F;

  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {

    a_posteriori_snr = 10.F * log10f(spectrum[k] / noise_spectrum[k]);

    if (a_posteriori_snr >= 0.F && a_posteriori_snr <= 20.F) {
      oversustraction_factor = -0.05F * (a_posteriori_snr) + 5.F;
    } else if (a_posteriori_snr < 0.F) {
      oversustraction_factor = 5.F;
    } else if (a_posteriori_snr > 20.F) {
      oversustraction_factor = 1.F;
    }

    noise_spectrum[k] *= parameters.noise_rescale * oversustraction_factor;
  }
}

static void masking_thresholds(OversubtractionCriterias *self,
                               const float *spectrum, float *noise_spectrum,
                               OversustractionParameters parameters) {

  if (parameters.masking_ceil > 1.F) {
    for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
      self->clean_signal_estimation[k] =
          fmaxf(spectrum[k] - noise_spectrum[k], FLT_MIN);
    }

    compute_masking_thresholds(self->masking_estimation, spectrum,
                               self->masking_thresholds);

    float max_masked_tmp =
        max_spectral_value(self->masking_thresholds, self->real_spectrum_size);
    float min_masked_tmp =
        min_spectral_value(self->masking_thresholds, self->real_spectrum_size);

    for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
      if (self->masking_thresholds[k] == max_masked_tmp) {
        self->alpha[k] = ALPHA_MIN;
        self->beta[k] = BETA_MIN;
      }
      if (self->masking_thresholds[k] == min_masked_tmp) {
        self->alpha[k] = parameters.masking_ceil;
        self->beta[k] = parameters.masking_floor;
      }
      if (self->masking_thresholds[k] < max_masked_tmp &&
          self->masking_thresholds[k] > min_masked_tmp) {
        const float normalized_value =
            (self->masking_thresholds[k] - min_masked_tmp) /
            (max_masked_tmp - min_masked_tmp);

        self->alpha[k] = (1.F - normalized_value) * ALPHA_MIN +
                         normalized_value * parameters.masking_ceil;
        self->beta[k] = (1.F - normalized_value) * BETA_MIN +
                        normalized_value * parameters.masking_floor;
      }
    }
  } else {
    initialize_spectrum_with_value(self->alpha, self->real_spectrum_size, 1.F);
  }

  for (uint32_t k = 1U; k < self->real_spectrum_size; k++) {
    noise_spectrum[k] *= parameters.noise_rescale * self->alpha[k];
  }
}