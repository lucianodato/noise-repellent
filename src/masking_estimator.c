/*
noise-repellent -- Noise Reduction LV2

Copyright 2016 Luciano Dato <lucianodato@gmail.com>

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
#include "spectral_utils.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N_BARK_BANDS 25
#define AT_SINE_WAVE_FREQ 1000.f
#define REFERENCE_LEVEL 90.f

#define BIAS 0
#define HIGH_FREQ_BIAS 20.f
#define S_AMP 1.f

#if BIAS
static const float relative_thresholds[N_BARK_BANDS] = {
    -16.f, -17.f, -18.f, -19.f, -20.f, -21.f, -22.f, -23.f, -24.f,
    -25.f, -25.f, -25.f, -25.f, -25.f, -25.f, -24.f, -23.f, -22.f,
    -19.f, -18.f, -18.f, -18.f, -18.f, -18.f, -18.f};
#endif

static void compute_spectral_spreading_function(MaskingEstimator *self);
static void convolve_with_spectral_spreading_function(MaskingEstimator *self);
static void compute_bark_mapping(MaskingEstimator *self);
static void compute_absolute_thresholds(MaskingEstimator *self);
static void spl_reference(MaskingEstimator *self);
static void compute_bark_spectrum(MaskingEstimator *self,
                                  const float *spectrum);
static float compute_tonality_factor(MaskingEstimator *self,
                                     const float *spectrum, uint32_t band);
static void convert_to_dbspl(MaskingEstimator *self, float *masking_thresholds);

struct MaskingEstimator {

  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;

  float *bark_z_spectrum;
  float *absolute_thresholds;
  float *spl_reference_values;
  float *input_fft_buffer_at;
  float *output_fft_buffer_at;

  float *spectral_spreading_function;
  float *unity_gain_bark_spectrum;
  float *spreaded_unity_gain_bark_spectrum;
  float *uint32_termediate_band_bins;
  float *n_bins_per_band;
  float *bark_spectrum;
  float *threshold_j;
  float *masking_offset;
  float *spreaded_spectrum;

  float *sinewave;
  float *window;
  float *fft_power_at;
  float *fft_power_at_dbspl;

  fftwf_plan forward_fft;
};

MaskingEstimator *masking_estimation_initialize(const uint32_t fft_size,
                                                const uint32_t sample_rate) {

  MaskingEstimator *self =
      (MaskingEstimator *)calloc(1, sizeof(MaskingEstimator));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->sample_rate = sample_rate;

  self->absolute_thresholds =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->bark_z_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->spl_reference_values =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->input_fft_buffer_at = (float *)calloc((self->fft_size), sizeof(float));
  self->output_fft_buffer_at = (float *)calloc((self->fft_size), sizeof(float));

  self->spectral_spreading_function =
      (float *)calloc((N_BARK_BANDS * N_BARK_BANDS), sizeof(float));
  self->unity_gain_bark_spectrum = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->spreaded_unity_gain_bark_spectrum =
      (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->uint32_termediate_band_bins =
      (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->n_bins_per_band = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->bark_spectrum = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->threshold_j = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->masking_offset = (float *)calloc(N_BARK_BANDS, sizeof(float));
  self->spreaded_spectrum = (float *)calloc(N_BARK_BANDS, sizeof(float));

  self->sinewave = (float *)calloc(self->fft_size, sizeof(float));
  self->window = (float *)calloc(self->fft_size, sizeof(float));
  self->fft_power_at =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->fft_power_at_dbspl =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->forward_fft =
      fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer_at,
                        self->output_fft_buffer_at, FFTW_R2HC, FFTW_ESTIMATE);

  compute_bark_mapping(self);
  compute_absolute_thresholds(self);
  spl_reference(self);
  compute_spectral_spreading_function(self);
  convolve_with_spectral_spreading_function(self);

  return self;
}

void masking_estimation_free(MaskingEstimator *self) {
  free(self->absolute_thresholds);
  free(self->bark_z_spectrum);
  free(self->spl_reference_values);
  free(self->input_fft_buffer_at);
  free(self->output_fft_buffer_at);

  free(self->spectral_spreading_function);
  free(self->unity_gain_bark_spectrum);
  free(self->spreaded_unity_gain_bark_spectrum);
  free(self->uint32_termediate_band_bins);
  free(self->n_bins_per_band);
  free(self->bark_spectrum);
  free(self->threshold_j);
  free(self->masking_offset);
  free(self->spreaded_spectrum);

  free(self->sinewave);
  free(self->window);
  free(self->fft_power_at);
  free(self->fft_power_at_dbspl);

  fftwf_destroy_plan(self->forward_fft);
  free(self);
}

void compute_masking_thresholds(MaskingEstimator *self, const float *spectrum,
                                float *masking_thresholds) {
  compute_bark_spectrum(self, spectrum);

  convolve_with_spectral_spreading_function(self);

  for (uint32_t j = 0; j < N_BARK_BANDS; j++) {

    const float tonality_factor = compute_tonality_factor(self, spectrum, j);

    self->masking_offset[j] = (tonality_factor * (14.5f + (float)(j + 1)) +
                               5.5f * (1.f - tonality_factor));

#if BIAS

    masking_offset[j] = relative_thresholds[j];

    if (j > 15) {
      masking_offset[j] += HIGH_FREQ_BIAS;
    }
#endif

    self->threshold_j[j] = powf(10.f, log10f(self->spreaded_spectrum[j]) -
                                          (self->masking_offset[j] / 10.f));

    self->threshold_j[j] -=
        10.f * log10f(self->spreaded_unity_gain_bark_spectrum[j]);

    float start_pos = 0.f;
    if (j == 0) {
      start_pos = 0;
    } else {
      start_pos = self->uint32_termediate_band_bins[j - 1];
    }
    const float end_pos = self->uint32_termediate_band_bins[j];

    for (uint32_t k = start_pos; k < end_pos; k++) {
      masking_thresholds[k] = self->threshold_j[j];
    }
  }

  convert_to_dbspl(self, masking_thresholds);

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    masking_thresholds[k] =
        fmaxf(masking_thresholds[k], self->absolute_thresholds[k]);
  }
}

static float bin_to_freq(const uint32_t bin_index, const float sample_rate,
                         const uint32_t half_fft_size) {
  return (float)bin_index * (sample_rate / half_fft_size / 2.f);
}

static void compute_bark_mapping(MaskingEstimator *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    const float frequency = (float)self->sample_rate /
                            (2.f * (float)(self->half_fft_size) * (float)k);
    self->bark_z_spectrum[k] = 1.f + 13.f * atanf(0.00076f * frequency) +
                               3.5f * atanf(powf(frequency / 7500.f, 2.f));
  }
}

static void compute_absolute_thresholds(MaskingEstimator *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {

    const float frequency =
        bin_to_freq(k, self->sample_rate, self->half_fft_size);
    self->absolute_thresholds[k] =
        3.64f * powf((frequency / 1000.f), -0.8f) -
        6.5f * expf(-0.6f * powf((frequency / 1000.f - 3.3f), 2.f)) +
        powf(10.f, -3.f) * powf((frequency / 1000.f), 4.f);
  }
}

static void spl_reference(MaskingEstimator *self) {
  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->sinewave[k] = S_AMP * sinf((2.f * M_PI * k * AT_SINE_WAVE_FREQ) /
                                     (float)self->sample_rate);
  }

  get_fft_window(self->window, self->fft_size, HANN_WINDOW);

  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->input_fft_buffer_at[k] = self->sinewave[k] * self->window[k];
  }

  fftwf_execute(self->forward_fft);

  get_fft_power_spectrum(self->output_fft_buffer_at, self->fft_size,
                         self->fft_power_at, self->half_fft_size);

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->fft_power_at_dbspl[k] =
        REFERENCE_LEVEL - 10.f * log10f(self->fft_power_at[k]);
  }

  memcpy(self->spl_reference_values, self->fft_power_at_dbspl,
         sizeof(float) * (self->half_fft_size + 1));
}

static void compute_spectral_spreading_function(MaskingEstimator *self) {
  for (uint32_t i = 0; i < N_BARK_BANDS; i++) {
    for (uint32_t j = 0; j < N_BARK_BANDS; j++) {
      const float y = (i + 1) - (j + 1);

      self->spectral_spreading_function[i * N_BARK_BANDS + j] =
          15.81f + 7.5f * (y + 0.474f) -
          17.5f * sqrtf(1.f + (y + 0.474f) * (y + 0.474f));

      self->spectral_spreading_function[i * N_BARK_BANDS + j] = powf(
          10.f, self->spectral_spreading_function[i * N_BARK_BANDS + j] / 10.f);
    }
  }
}

static void convolve_with_spectral_spreading_function(MaskingEstimator *self) {
  for (uint32_t i = 0; i < N_BARK_BANDS; i++) {
    self->spreaded_spectrum[i] = 0.f;
    for (uint32_t j = 0; j < N_BARK_BANDS; j++) {
      self->spreaded_spectrum[i] +=
          self->spectral_spreading_function[i * N_BARK_BANDS + j] *
          self->bark_spectrum[j];
    }
  }
}

static void compute_bark_spectrum(MaskingEstimator *self,
                                  const float *spectrum) {
  uint32_t last_position = 0;

  for (uint32_t j = 0; j < N_BARK_BANDS; j++) {
    uint32_t counter = 0;
    if (j == 0) {
      counter = 1;
    }

    self->bark_spectrum[j] = 0.f;

    while (floorf(self->bark_z_spectrum[last_position + counter]) == (j + 1)) {
      self->bark_spectrum[j] += spectrum[last_position + counter];
      counter++;
    }

    last_position += counter;

    self->n_bins_per_band[j] = counter;
    self->uint32_termediate_band_bins[j] = last_position;
  }
}

static void convert_to_dbspl(MaskingEstimator *self,
                             float *masking_thresholds) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    masking_thresholds[k] += self->spl_reference_values[k];
  }
}

static float compute_tonality_factor(MaskingEstimator *self,
                                     const float *spectrum, uint32_t band) {
  float sum_bins = 0.f;
  float sum_log_bins = 0.f;
  uint32_t start_pos = 0;
  uint32_t end_pos = 0;

  if (band == 0) {
    start_pos = band;
    end_pos = self->n_bins_per_band[band];
  } else {
    start_pos = self->uint32_termediate_band_bins[band - 1];
    end_pos = self->uint32_termediate_band_bins[band - 1] +
              self->n_bins_per_band[band];
  }

  for (uint32_t k = start_pos; k < end_pos; k++) {
    sum_bins += spectrum[k];
    sum_log_bins += log10f(spectrum[k]);
  }

  const float SFM = 10.f * (sum_log_bins / (self->n_bins_per_band[band]) -
                            log10f(sum_bins / (self->n_bins_per_band[band])));

  const float tonality_factor = fminf(SFM / -60.f, 1.f);

  return tonality_factor;
}
