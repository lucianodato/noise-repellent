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
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define WHITENING_DECAY_RATE 1000.f
#define WHITENING_FLOOR 0.02f

static void fft_denoiser_update_wetdry_target(SpectralDenoiser *self);
static void fft_denoiser_soft_bypass(SpectralDenoiser *self);
static void get_denoised_spectrum(SpectralDenoiser *self);
static void get_residual_spectrum(SpectralDenoiser *self);
static void get_final_spectrum(SpectralDenoiser *self);

typedef struct {
  float tau;
  float wet_dry_target;
  float wet_dry;
} SoftBypass;

typedef struct {
  float *residual_max_spectrum;
  float *whitened_residual_spectrum;
  float max_decay_rate;
  uint32_t whitening_window_count;
} Whitening;

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

  SoftBypass crossfade_spectrum;
  Whitening whiten_spectrum;
  SpectralDenoiseBuilder denoise_builder;

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

  self->crossfade_spectrum.tau =
      (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->sample_rate));
  self->crossfade_spectrum.wet_dry = 0.f;

  self->whiten_spectrum.whitened_residual_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->whiten_spectrum.residual_max_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->whiten_spectrum.max_decay_rate =
      expf(-1000.f / (((WHITENING_DECAY_RATE)*self->sample_rate) / self->hop));
  self->whiten_spectrum.whitening_window_count = 0.f;

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

  free(self->fft_spectrum);
  free(self->processed_fft_spectrum);
  free(self->denoise_builder.gain_spectrum);
  free(self->denoise_builder.residual_spectrum);
  free(self->denoise_builder.denoised_spectrum);
  free(self->whiten_spectrum.whitened_residual_spectrum);
  free(self->whiten_spectrum.residual_max_spectrum);
  free(self->processing_spectrums.power_spectrum);
  free(self);
}

void spectral_denoiser_run(SPECTRAL_PROCESSOR instance, float *fft_spectrum) {
  SpectralDenoiser *self = (SpectralDenoiser *)instance;

  fft_denoiser_update_wetdry_target(self);

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

  fft_denoiser_soft_bypass(self);

  memcpy(fft_spectrum, self->processed_fft_spectrum,
         sizeof(float) * self->half_fft_size + 1);
}

static void residual_spectrum_whitening(SpectralDenoiser *self,
                                        const float whitening_factor) {
  self->whiten_spectrum.whitening_window_count++;

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->whiten_spectrum.whitening_window_count > 1.f) {
      self->whiten_spectrum.residual_max_spectrum[k] = fmaxf(
          fmaxf(self->denoise_builder.residual_spectrum[k], WHITENING_FLOOR),
          self->whiten_spectrum.residual_max_spectrum[k] *
              self->whiten_spectrum.max_decay_rate);
    } else {
      self->whiten_spectrum.residual_max_spectrum[k] =
          fmaxf(self->denoise_builder.residual_spectrum[k], WHITENING_FLOOR);
    }
  }

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->denoise_builder.residual_spectrum[k] > FLT_MIN) {
      self->whiten_spectrum.whitened_residual_spectrum[k] =
          self->denoise_builder.residual_spectrum[k] /
          self->whiten_spectrum.residual_max_spectrum[k];

      self->denoise_builder.residual_spectrum[k] =
          (1.f - whitening_factor) *
              self->denoise_builder.residual_spectrum[k] +
          whitening_factor *
              self->whiten_spectrum.whitened_residual_spectrum[k];
    }
  }
}

static void fft_denoiser_update_wetdry_target(SpectralDenoiser *self) {
  if (self->denoise_parameters->enable) {
    self->crossfade_spectrum.wet_dry_target = 1.f;
  } else {
    self->crossfade_spectrum.wet_dry_target = 0.f;
  }

  self->crossfade_spectrum.wet_dry +=
      self->crossfade_spectrum.tau * (self->crossfade_spectrum.wet_dry_target -
                                      self->crossfade_spectrum.wet_dry) +
      FLT_MIN;
}

static void fft_denoiser_soft_bypass(SpectralDenoiser *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->processed_fft_spectrum[k] =
        (1.f - self->crossfade_spectrum.wet_dry) * self->fft_spectrum[k] +
        self->processed_fft_spectrum[k] * self->crossfade_spectrum.wet_dry;
  }
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
    residual_spectrum_whitening(self,
                                self->denoise_parameters->whitening_factor);
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
