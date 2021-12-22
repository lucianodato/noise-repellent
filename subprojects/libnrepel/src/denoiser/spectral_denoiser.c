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

#include "spectral_denoiser.h"
#include "../shared/spectral_utils.h"
#include "gain_estimator.h"
#include "spectral_whitening.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void denoise_build(SpectralDenoiser *self, float *fft_spectrum);

struct SpectralDenoiser {
  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;
  uint32_t hop;

  float *gain_spectrum;
  float *residual_spectrum;
  float *denoised_spectrum;

  SpectralWhitening *whitener;
  NoiseProfile *noise_profile;
  GainEstimator *gain_estimation;
  ProcessorParameters *denoise_parameters;
};

SpectralDenoiserHandle spectral_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile *noise_profile,
    ProcessorParameters *parameters) {
  SpectralDenoiser *self =
      (SpectralDenoiser *)calloc(1U, sizeof(SpectralDenoiser));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2U;
  self->hop = self->fft_size / overlap_factor;
  self->sample_rate = sample_rate;

  self->gain_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));
  initialize_spectrum_to_ones(self->gain_spectrum, self->half_fft_size + 1U);

  self->noise_profile = noise_profile;
  self->denoise_parameters = parameters;

  self->gain_estimation = gain_estimation_initialize(
      self->fft_size, self->sample_rate, self->hop, self->denoise_parameters);

  self->residual_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));
  self->denoised_spectrum =
      (float *)calloc((self->half_fft_size + 1U), sizeof(float));

  self->whitener = spectral_whitening_initialize(self->fft_size,
                                                 self->sample_rate, self->hop);

  return self;
}

void spectral_denoiser_free(SpectralDenoiserHandle instance) {
  SpectralDenoiser *self = (SpectralDenoiser *)instance;

  gain_estimation_free(self->gain_estimation);
  spectral_whitening_free(self->whitener);

  free(self->residual_spectrum);
  free(self->denoised_spectrum);
  free(self->gain_spectrum);
  free(self);
}

bool spectral_denoiser_run(SpectralDenoiserHandle instance,
                           float *fft_spectrum) {
  if (!fft_spectrum || !instance) {
    return false;
  }

  SpectralDenoiser *self = (SpectralDenoiser *)instance;

  gain_estimation_run(self->gain_estimation, fft_spectrum,
                      get_noise_profile(self->noise_profile),
                      self->gain_spectrum);

  denoise_build(self, fft_spectrum);

  return true;
}

static void get_denoised_spectrum(SpectralDenoiser *self,
                                  const float *fft_spectrum) {
  for (uint32_t k = 1U; k <= self->half_fft_size; k++) {
    self->denoised_spectrum[k] = fft_spectrum[k] * self->gain_spectrum[k];
  }
}

static void get_residual_spectrum(SpectralDenoiser *self,
                                  const float *fft_spectrum) {
  for (uint32_t k = 1U; k <= self->half_fft_size; k++) {
    self->residual_spectrum[k] = fft_spectrum[k] - self->denoised_spectrum[k];
  }

  if (self->denoise_parameters->whitening_factor > 0.f) {
    spectral_whitening_run(self->whitener,
                           self->denoise_parameters->whitening_factor,
                           self->residual_spectrum);
  }
}

static void denoise_build(SpectralDenoiser *self, float *fft_spectrum) {

  get_denoised_spectrum(self, fft_spectrum);

  get_residual_spectrum(self, fft_spectrum);

  if (self->denoise_parameters->residual_listen) {
    for (uint32_t k = 1U; k <= self->half_fft_size; k++) {
      fft_spectrum[k] = self->residual_spectrum[k];
    }
  } else {
    for (uint32_t k = 1U; k <= self->half_fft_size; k++) {
      fft_spectrum[k] = self->denoised_spectrum[k] +
                        self->residual_spectrum[k] *
                            self->denoise_parameters->reduction_amount;
    }
  }
}