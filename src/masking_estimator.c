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
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

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

void compute_spectral_spreading_function(MaskingEstimator *self);
void convolve_with_spectral_spreading_function(MaskingEstimator *self,
                                               const float *bark_spectrum,
                                               float *spreaded_spectrum);
void compute_bark_mapping(MaskingEstimator *self);
void compute_absolute_thresholds(MaskingEstimator *self);
void spl_reference(MaskingEstimator *self);

struct MaskingEstimator {

  int fft_size;
  int half_fft_size;
  int samp_rate;

  float *bark_z_spectrum;
  float *absolute_thresholds;
  float *spl_reference_values;
  float *input_fft_buffer_at;
  float *output_fft_buffer_at;
  float *spectral_spreading_function;
  float *unity_gain_bark_spectrum;
  float *spreaded_unity_gain_bark_spectrum;
  fftwf_plan forward_fft;
};

MaskingEstimator *masking_estimation_initialize(int fft_size, int samp_rate) {

  MaskingEstimator *self =
      (MaskingEstimator *)calloc(1, sizeof(MaskingEstimator));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->samp_rate = samp_rate;

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

  self->forward_fft =
      fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer_at,
                        self->output_fft_buffer_at, FFTW_R2HC, FFTW_ESTIMATE);

  compute_bark_mapping(self);
  compute_absolute_thresholds(self);
  spl_reference(self);
  compute_spectral_spreading_function(self);
  convolve_with_spectral_spreading_function(
      self, self->unity_gain_bark_spectrum,
      self->spreaded_unity_gain_bark_spectrum);

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
  fftwf_destroy_plan(self->forward_fft);
  free(self);
}

float bin_to_freq(int i, float samp_rate, int N) {
  return (float)i * (samp_rate / N / 2.f);
}

void compute_bark_mapping(MaskingEstimator *self) {
  for (int k = 1; k <= self->half_fft_size; k++) {
    float freq = 0.f;
    freq = (float)self->samp_rate /
           (2.f * (float)(self->half_fft_size) * (float)k);
    self->bark_z_spectrum[k] = 1.f + 13.f * atanf(0.00076f * freq) +
                               3.5f * atanf(powf(freq / 7500.f, 2.f));
  }
}

void compute_absolute_thresholds(MaskingEstimator *self) {
  for (int k = 1; k <= self->half_fft_size; k++) {
    float freq = 0.f;

    freq = bin_to_freq(k, self->samp_rate, self->half_fft_size);
    self->absolute_thresholds[k] =
        3.64f * powf((freq / 1000.f), -0.8f) -
        6.5f * expf(-0.6f * powf((freq / 1000.f - 3.3f), 2.f)) +
        powf(10.f, -3.f) * powf((freq / 1000.f), 4.f);
  }
}

void hanning_window(float *window, int N) {
  for (int k = 0; k < N; k++) {
    float p = ((float)(k)) / ((float)(N));
    window[k] = 0.5 - 0.5 * cosf(2.f * M_PI * p);
  }
}

void get_power_spectrum(MaskingEstimator *self, float *window,
                        const float *signal, float *power_spectrum) {
  hanning_window(window, self->fft_size);
  for (int k = 0; k < self->fft_size; k++) {
    self->input_fft_buffer_at[k] = signal[k] * window[k];
  }

  fftwf_execute(self->forward_fft);

  float real_bin = self->output_fft_buffer_at[0];

  power_spectrum[0] = real_bin * real_bin;

  for (int k = 1; k <= self->half_fft_size; k++) {
    float imag_bin = 0.f;
    float p2 = 0.f;

    real_bin = self->output_fft_buffer_at[k];
    imag_bin = self->output_fft_buffer_at[self->fft_size - k];

    if (k < self->half_fft_size) {
      p2 = (real_bin * real_bin + imag_bin * imag_bin);
    } else {

      p2 = real_bin * real_bin;
    }

    power_spectrum[k] = p2;
  }
}

void spl_reference(MaskingEstimator *self) {
  float sinewave[self->fft_size];
  float window[self->fft_size];
  float fft_power_at[self->half_fft_size + 1];
  float fft_power_at_dbspl[self->half_fft_size + 1];

  for (int k = 0; k < self->fft_size; k++) {
    sinewave[k] = S_AMP * sinf((2.f * M_PI * k * AT_SINE_WAVE_FREQ) /
                               (float)self->samp_rate);
  }

  get_power_spectrum(self, window, sinewave, fft_power_at);

  for (int k = 1; k <= self->half_fft_size; k++) {
    fft_power_at_dbspl[k] = REFERENCE_LEVEL - 10.f * log10f(fft_power_at[k]);
  }

  memcpy(self->spl_reference_values, fft_power_at_dbspl,
         sizeof(float) * (self->half_fft_size + 1));
}

void compute_spectral_spreading_function(MaskingEstimator *self) {
  float y = 0.f;
  for (int i = 0; i < N_BARK_BANDS; i++) {
    for (int j = 0; j < N_BARK_BANDS; j++) {
      y = (i + 1) - (j + 1);

      self->spectral_spreading_function[i * N_BARK_BANDS + j] =
          15.81f + 7.5f * (y + 0.474f) -
          17.5f * sqrtf(1.f + (y + 0.474f) * (y + 0.474f));

      self->spectral_spreading_function[i * N_BARK_BANDS + j] = powf(
          10.f, self->spectral_spreading_function[i * N_BARK_BANDS + j] / 10.f);
    }
  }
}

void convolve_with_spectral_spreading_function(MaskingEstimator *self,
                                               const float *bark_spectrum,
                                               float *spreaded_spectrum) {
  for (int i = 0; i < N_BARK_BANDS; i++) {
    spreaded_spectrum[i] = 0.f;
    for (int j = 0; j < N_BARK_BANDS; j++) {
      spreaded_spectrum[i] +=
          self->spectral_spreading_function[i * N_BARK_BANDS + j] *
          bark_spectrum[j];
    }
  }
}

void compute_bark_spectrum(MaskingEstimator *self, float *bark_spectrum,
                           const float *spectrum, float *intermediate_band_bins,
                           float *n_bins_per_band) {
  int last_position = 0;

  for (int j = 0; j < N_BARK_BANDS; j++) {
    int cont = 0;
    if (j == 0) {
      cont = 1;
    }

    bark_spectrum[j] = 0.f;

    while (floorf(self->bark_z_spectrum[last_position + cont]) == (j + 1)) {
      bark_spectrum[j] += spectrum[last_position + cont];
      cont++;
    }

    last_position += cont;

    n_bins_per_band[j] = cont;
    intermediate_band_bins[j] = last_position;
  }
}

void convert_to_dbspl(MaskingEstimator *self, float *masking_thresholds) {
  for (int k = 1; k <= self->half_fft_size; k++) {
    masking_thresholds[k] += self->spl_reference_values[k];
  }
}

float compute_tonality_factor(float *spectrum,
                              const float *intermediate_band_bins,
                              const float *n_bins_per_band, int band) {
  int k = 0;
  float SFM = 0.f;
  float tonality_factor = 0.f;
  float sum_bins = 0.f;
  float sum_log_bins = 0.f;
  int start_pos = 0;
  int end_pos = 0;

  if (band == 0) {
    start_pos = band;
    end_pos = n_bins_per_band[band];
  } else {
    start_pos = intermediate_band_bins[band - 1];
    end_pos = intermediate_band_bins[band - 1] + n_bins_per_band[band];
  }

  for (k = start_pos; k < end_pos; k++) {
    sum_bins += spectrum[k];
    sum_log_bins += log10f(spectrum[k]);
  }

  SFM = 10.f * (sum_log_bins / (n_bins_per_band[band]) -
                log10f(sum_bins / (n_bins_per_band[band])));

  tonality_factor = fminf(SFM / -60.f, 1.f);

  return tonality_factor;
}

void compute_masking_thresholds(MaskingEstimator *self, float *spectrum,
                                float *masking_thresholds) {
  float intermediate_band_bins[N_BARK_BANDS];
  float n_bins_per_band[N_BARK_BANDS];
  float bark_spectrum[N_BARK_BANDS];
  float threshold_j[N_BARK_BANDS];
  float masking_offset[N_BARK_BANDS];
  float spreaded_spectrum[N_BARK_BANDS];

  compute_bark_spectrum(self, bark_spectrum, spectrum, intermediate_band_bins,
                        n_bins_per_band);

  convolve_with_spectral_spreading_function(self, bark_spectrum,
                                            spreaded_spectrum);

  for (int j = 0; j < N_BARK_BANDS; j++) {
    float tonality_factor = 0.f;
    float start_pos = 0.f;
    float end_pos = 0.f;

    tonality_factor = compute_tonality_factor(spectrum, intermediate_band_bins,
                                              n_bins_per_band, j);

    masking_offset[j] = (tonality_factor * (14.5f + (float)(j + 1)) +
                         5.5f * (1.f - tonality_factor));

#if BIAS

    masking_offset[j] = relative_thresholds[j];

    if (j > 15) {
      masking_offset[j] += HIGH_FREQ_BIAS;
    }
#endif

    threshold_j[j] =
        powf(10.f, log10f(spreaded_spectrum[j]) - (masking_offset[j] / 10.f));

    threshold_j[j] -= 10.f * log10f(self->spreaded_unity_gain_bark_spectrum[j]);

    if (j == 0) {
      start_pos = 0;
    } else {
      start_pos = intermediate_band_bins[j - 1];
    }
    end_pos = intermediate_band_bins[j];

    for (int k = start_pos; k < end_pos; k++) {
      masking_thresholds[k] = threshold_j[j];
    }
  }

  convert_to_dbspl(self, masking_thresholds);

  for (int k = 1; k <= self->half_fft_size; k++) {
    masking_thresholds[k] =
        fmaxf(masking_thresholds[k], self->absolute_thresholds[k]);
  }
}