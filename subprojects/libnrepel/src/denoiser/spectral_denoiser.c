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

#include "spectral_denoiser.h"
#include "../shared/spectral_features.h"
#include "gain_estimator.h"
#include "spectral_whitening.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void get_denoised_spectrum(SpectralDenoiser *self);
static void get_residual_spectrum(SpectralDenoiser *self);
static void get_final_spectrum(SpectralDenoiser *self);

typedef struct {
  float *gain_spectrum;
  float *residual_spectrum;
  float *denoised_spectrum;
} SpectralDenoiseBuilder;

struct SpectralDenoiser {
  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;
  uint32_t hop;

  float *fft_spectrum;
  float *processed_fft_spectrum;

  NoiseProfile *noise_profile;

  SpectralDenoiseBuilder denoise_builder;

  SpectralWhitening *whitener;
  GainEstimator *gain_estimation;
  ProcessorParameters *denoise_parameters;

  SpectralFeatures processing_spectrums;
};

SpectralDenoiser *spectral_denoiser_initialize(
    const uint32_t sample_rate, const uint32_t fft_size,
    const uint32_t overlap_factor, NoiseProfile *noise_profile,
    ProcessorParameters *parameters) {
  SpectralDenoiser *self =
      (SpectralDenoiser *)calloc(1, sizeof(SpectralDenoiser));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->hop = self->fft_size / overlap_factor;
  self->sample_rate = sample_rate;

  self->fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));
  self->processed_fft_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->denoise_builder.residual_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->denoise_builder.denoised_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->denoise_builder.gain_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->whitener = spectral_whitening_initialize(self->fft_size,
                                                 self->sample_rate, self->hop);

  self->noise_profile = noise_profile;
  self->denoise_parameters = parameters;

  self->gain_estimation = gain_estimation_initialize(
      self->fft_size, self->sample_rate, self->hop, self->denoise_parameters);

  self->processing_spectrums.power_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  return self;
}

void spectral_denoiser_free(SpectralDenoiser *self) {
  gain_estimation_free(self->gain_estimation);
  spectral_whitening_free(self->whitener);

  free(self->fft_spectrum);
  free(self->processed_fft_spectrum);
  free(self->denoise_builder.gain_spectrum);
  free(self->denoise_builder.residual_spectrum);
  free(self->denoise_builder.denoised_spectrum);
  free(self->processing_spectrums.power_spectrum);
  free(self);
}

void spectral_denoiser_run(SPECTRAL_PROCESSOR instance, float *fft_spectrum) {
  SpectralDenoiser *self = (SpectralDenoiser *)instance;

  memcpy(self->fft_spectrum, fft_spectrum,
         sizeof(float) * self->half_fft_size + 1);

  get_fft_power_spectrum(self->fft_spectrum, self->fft_size,
                         self->processing_spectrums.power_spectrum,
                         self->half_fft_size);

  gain_estimation_run(self->gain_estimation,
                      self->processing_spectrums.power_spectrum,
                      get_noise_profile(self->noise_profile),
                      self->denoise_builder.gain_spectrum);

  get_denoised_spectrum(self);

  get_residual_spectrum(self);

  get_final_spectrum(self);

  memcpy(fft_spectrum, self->processed_fft_spectrum,
         sizeof(float) * self->half_fft_size + 1);
}

static void get_denoised_spectrum(SpectralDenoiser *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->denoise_builder.denoised_spectrum[k] =
        self->fft_spectrum[k] * self->denoise_builder.gain_spectrum[k];
  }
}

static void get_residual_spectrum(SpectralDenoiser *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->denoise_builder.residual_spectrum[k] =
        self->fft_spectrum[k] - self->denoise_builder.denoised_spectrum[k];
  }

  if (self->denoise_parameters->whitening_factor > 0.f) {
    spectral_whitening_run(self->whitener,
                           self->denoise_parameters->whitening_factor,
                           self->denoise_builder.residual_spectrum);
  }
}

static void get_final_spectrum(SpectralDenoiser *self) {
  if (self->denoise_parameters->residual_listen) {
    for (uint32_t k = 1; k <= self->half_fft_size; k++) {
      self->processed_fft_spectrum[k] =
          self->denoise_builder.residual_spectrum[k];
    }
  } else {
    for (uint32_t k = 1; k <= self->half_fft_size; k++) {
      self->processed_fft_spectrum[k] =
          self->denoise_builder.denoised_spectrum[k] +
          self->denoise_builder.residual_spectrum[k] *
              self->denoise_parameters->reduction_amount;
    }
  }
}
