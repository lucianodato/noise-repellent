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

#include "postfilter.h"
#include "configurations.h"
#include "fft_transform.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct PostFilter {
  float *postfilter;
  FftTransform *gain_fft_spectrum;
  FftTransform *postfilter_fft_spectrum;
  bool preserve_minimun;

  float *pf_gain_spectrum;

  float snr_threshold;
  float postfilter_scale;
  uint32_t fft_size;
  uint32_t real_spectrum_size;
};

PostFilter *postfilter_initialize(const uint32_t fft_size) {
  PostFilter *self = (PostFilter *)calloc(1U, sizeof(PostFilter));

  self->fft_size = fft_size;
  self->real_spectrum_size = self->fft_size / 2U + 1U;
  self->snr_threshold = POSTFILTER_THRESHOLD;
  self->postfilter_scale = POSTFILTER_SCALE;
  self->preserve_minimun = (bool)PRESERVE_MINIMUN_GAIN;

  self->gain_fft_spectrum = fft_transform_initialize_bins(self->fft_size);
  self->postfilter_fft_spectrum = fft_transform_initialize_bins(self->fft_size);
  self->pf_gain_spectrum = (float *)calloc(self->fft_size, sizeof(float));

  self->postfilter = (float *)calloc(self->fft_size, sizeof(float));

  return self;
}

void postfilter_free(PostFilter *self) {
  fft_transform_free(self->gain_fft_spectrum);
  fft_transform_free(self->postfilter_fft_spectrum);

  free(self->postfilter);
  free(self->pf_gain_spectrum);

  free(self);
}

static void calculate_postfilter(PostFilter *self, const float *spectrum,
                                 const float *gain_spectrum) {
  float clean_signal_sum = 0.F;
  float noisy_signa_sum = 0.F;
  float a_priori_snr = 0.F;
  float lambda = 0.F;
  float threshold_decision = 0.F;

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    clean_signal_sum += powf(spectrum[k] * gain_spectrum[k], 2.F);
    noisy_signa_sum += powf(spectrum[k], 2.F);
  }

  a_priori_snr = clean_signal_sum / noisy_signa_sum;

  if (a_priori_snr >= self->snr_threshold) {
    threshold_decision = 1.F;
  } else {
    threshold_decision = a_priori_snr;
  }

  if (threshold_decision == 1.F) {
    lambda = 1.F;
  } else {
    lambda = 2.F * roundf(self->postfilter_scale *
                          (1.F - threshold_decision / self->snr_threshold)) +
             1.F;
  }

  for (uint32_t k = 0U; k < self->real_spectrum_size; k++) {
    if ((float)k < lambda) {
      self->postfilter[k] = 1.F / lambda;
    } else {
      self->postfilter[k] = 0.F;
    }
  }
}

bool postfilter_apply(PostFilter *self, const float *spectrum,
                      float *gain_spectrum) {
  if (!spectrum || !gain_spectrum) {
    return false;
  }

  memcpy(self->pf_gain_spectrum, gain_spectrum, self->fft_size * sizeof(float));

  calculate_postfilter(self, spectrum, self->pf_gain_spectrum);

  fft_load_input_samples(self->gain_fft_spectrum, self->pf_gain_spectrum);
  fft_load_input_samples(self->postfilter_fft_spectrum, self->postfilter);

  compute_forward_fft(self->gain_fft_spectrum);
  compute_forward_fft(self->postfilter_fft_spectrum);

  for (uint32_t k = 0U; k < self->fft_size; k++) {
    get_fft_output_buffer(self->gain_fft_spectrum)[k] *=
        get_fft_output_buffer(self->postfilter_fft_spectrum)[k];
  }

  compute_backward_fft(self->gain_fft_spectrum);

  for (uint32_t k = 0U; k < self->fft_size; k++) {
    self->pf_gain_spectrum[k] =
        get_fft_input_buffer(self->gain_fft_spectrum)[k] /
        (float)self->fft_size;
  }

  if (self->preserve_minimun) {
    min_spectrum(gain_spectrum, self->pf_gain_spectrum, self->fft_size);
  } else {
    memcpy(gain_spectrum, self->pf_gain_spectrum,
           self->fft_size * sizeof(float));
  }

  return true;
}
