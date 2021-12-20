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
#include "denoise_builder.h"
#include "gain_estimator.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct SpectralDenoiser {
  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;
  uint32_t hop;

  float *gain_spectrum;

  NoiseProfile *noise_profile;
  DenoiseBuilder *denoise_builder;
  GainEstimator *gain_estimation;
  ProcessorParameters *denoise_parameters;
  SpectralFeatures *spectral_features;
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

  self->gain_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->denoise_builder = denoise_builder_initialize(
      self->half_fft_size, self->fft_size, self->hop, self->sample_rate);

  self->noise_profile = noise_profile;
  self->denoise_parameters = parameters;

  self->gain_estimation = gain_estimation_initialize(
      self->fft_size, self->sample_rate, self->hop, self->denoise_parameters);

  self->spectral_features = spectral_features_initialize(self->half_fft_size);

  return self;
}

void spectral_denoiser_free(SpectralDenoiser *self) {
  gain_estimation_free(self->gain_estimation);
  spectral_features_free(self->spectral_features);
  denoise_builder_free(self->denoise_builder);

  free(self->gain_spectrum);
  free(self);
}

void spectral_denoiser_run(SPECTRAL_PROCESSOR instance, float *fft_spectrum) {
  SpectralDenoiser *self = (SpectralDenoiser *)instance;

  compute_power_spectrum(
      self->spectral_features, fft_spectrum,
      self->fft_size); // TODO Move this inside gain estimation

  gain_estimation_run(
      self->gain_estimation, get_power_spectrum(self->spectral_features),
      get_noise_profile(self->noise_profile), self->gain_spectrum);

  denoise_build(self->denoise_builder,
                self->denoise_parameters->residual_listen,
                self->denoise_parameters->reduction_amount,
                self->denoise_parameters->whitening_factor, self->gain_spectrum,
                fft_spectrum);
}
