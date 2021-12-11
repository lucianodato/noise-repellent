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

#include "fft_denoiser.h"
#include "gain_estimator.h"
#include "noise_estimator.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define WHITENING_DECAY_RATE 1000.f
#define WHITENING_FLOOR 0.02f

static void get_info_from_bins(float *fft_power, float *fft_magnitude,
                               float *fft_phase, uint32_t half_fft_size,
                               uint32_t fft_size, const float *fft_buffer);
static bool is_empty(const float *spectrum, uint32_t half_fft_size);
static void fft_denoiser_update_wetdry_target(FFTDenoiser *self, bool enable);
static void fft_denoiser_soft_bypass(FFTDenoiser *self);
static void get_denoised_spectrum(FFTDenoiser *self);
static void get_residual_spectrum(FFTDenoiser *self, float whitening_factor);

static void get_final_spectrum(FFTDenoiser *self, bool residual_listen,
                               float reduction_amount);
static inline float from_db_to_coefficient(float gain_db);

struct FFTDenoiser {
  uint32_t fft_size;
  uint32_t half_fft_size;
  uint32_t sample_rate;
  uint32_t hop;

  float *fft_spectrum;
  float *processed_fft_spectrum;

  float tau;
  float wet_dry_target;
  float wet_dry;

  float *gain_spectrum;
  float *residual_spectrum;
  float *denoised_spectrum;
  float *whitened_residual_spectrum;

  float *power_spectrum;
  float *phase_spectrum;
  float *magnitude_spectrum;

  GainEstimator *gain_estimation;
  NoiseEstimator *noise_estimation;
  NoiseProfile *noise_profile;
  DenoiseParameters *denoise_parameters;

  float *residual_max_spectrum;
  float max_decay_rate;
  uint32_t whitening_window_count;
};

FFTDenoiser *fft_denoiser_initialize(const uint32_t sample_rate,
                                     const uint32_t fft_size,
                                     const uint32_t overlap_factor) {
  FFTDenoiser *self = (FFTDenoiser *)calloc(1, sizeof(FFTDenoiser));

  self->fft_size = fft_size;
  self->half_fft_size = self->fft_size / 2;
  self->hop = self->fft_size / overlap_factor;
  self->sample_rate = sample_rate;

  self->fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));
  self->processed_fft_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->sample_rate));
  self->wet_dry = 0.f;

  self->residual_max_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->max_decay_rate =
      expf(-1000.f / (((WHITENING_DECAY_RATE)*self->sample_rate) / self->hop));

  self->residual_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->denoised_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->gain_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->whitened_residual_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->whitening_window_count = 0.f;

  self->gain_estimation =
      gain_estimation_initialize(self->fft_size, self->sample_rate, self->hop);
  self->noise_estimation = noise_estimation_initialize(self->fft_size);

  self->power_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->magnitude_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));
  self->phase_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  return self;
}

void fft_denoiser_free(FFTDenoiser *self) {
  gain_estimation_free(self->gain_estimation);
  noise_estimation_free(self->noise_estimation);

  free(self->fft_spectrum);
  free(self->processed_fft_spectrum);
  free(self->gain_spectrum);
  free(self->residual_spectrum);
  free(self->whitened_residual_spectrum);
  free(self->denoised_spectrum);
  free(self->residual_max_spectrum);
  free(self->power_spectrum);
  free(self->magnitude_spectrum);
  free(self->phase_spectrum);
  free(self);
}

void load_denoise_parameters(FFTDenoiser *self,
                             DenoiseParameters *new_parameters) {
  self->denoise_parameters = new_parameters;
}

void load_noise_profile(FFTDenoiser *self, NoiseProfile *noise_profile) {
  self->noise_profile = noise_profile;
}

void fft_denoiser_run(FFTDenoiser *self, float *fft_spectrum) {

  const bool enable = (bool)*self->denoise_parameters->enable;
  const bool learn_noise = (bool)*self->denoise_parameters->learn_noise;
  const bool residual_listen = (bool)*self->denoise_parameters->residual_listen;
  const float transient_protection =
      *self->denoise_parameters->transient_threshold;
  const float masking =
      *self->denoise_parameters->masking_ceiling_limit / 100.f;
  const float release = *self->denoise_parameters->release_time;
  const float noise_rescale = *self->denoise_parameters->noise_rescale;
  const float reduction_amount = from_db_to_coefficient(
      *self->denoise_parameters->reduction_amount * -1.f);
  const float whitening_factor = *self->denoise_parameters->whitening_factor;
  float *noise_spectrum = self->noise_profile->noise_profile;

  fft_denoiser_update_wetdry_target(self, enable);

  memcpy(self->fft_spectrum, fft_spectrum, sizeof(float) * self->fft_size);

  get_info_from_bins(self->power_spectrum, self->magnitude_spectrum,
                     self->phase_spectrum, self->half_fft_size, self->fft_size,
                     self->fft_spectrum);

  if (!is_empty(self->power_spectrum, self->half_fft_size)) {
    if (learn_noise) {
      noise_estimation_run(self->noise_estimation, noise_spectrum,
                           self->power_spectrum);
    } else {
      if (is_noise_estimation_available(self->noise_estimation)) {
        gain_estimation_run(self->gain_estimation, self->power_spectrum,
                            noise_spectrum, self->gain_spectrum,
                            transient_protection, masking, release,
                            noise_rescale);

        get_denoised_spectrum(self);

        get_residual_spectrum(self, whitening_factor);

        get_final_spectrum(self, residual_listen, reduction_amount);
      }
    }
  }

  fft_denoiser_soft_bypass(self);

  memcpy(fft_spectrum, self->processed_fft_spectrum,
         sizeof(float) * self->half_fft_size + 1);
}

static void get_info_from_bins(float *fft_power, float *fft_magnitude,
                               float *fft_phase, const uint32_t half_fft_size,
                               const uint32_t fft_size,
                               const float *fft_buffer) {
  float real_bin = fft_buffer[0];

  fft_power[0] = real_bin * real_bin;
  fft_magnitude[0] = real_bin;
  fft_phase[0] = atan2f(real_bin, 0.f);

  for (uint32_t k = 1; k <= half_fft_size; k++) {
    float magnitude = 0.f;
    float power = 0.f;
    float phase = 0.f;

    real_bin = fft_buffer[k];
    float imag_bin = fft_buffer[fft_size - k];

    if (k < half_fft_size) {
      power = (real_bin * real_bin + imag_bin * imag_bin);
      magnitude = sqrtf(power);
      phase = atan2f(real_bin, imag_bin);
    } else {
      power = real_bin * real_bin;
      magnitude = real_bin;
      phase = atan2f(real_bin, 0.f);
    }

    fft_power[k] = power;
    fft_magnitude[k] = magnitude;
    fft_phase[k] = phase;
  }
}

static bool is_empty(const float *spectrum, const uint32_t half_fft_size) {
  for (uint32_t k = 1; k <= half_fft_size; k++) {
    if (spectrum[k] > FLT_MIN) {
      return false;
    }
  }
  return true;
}

static void fft_denoiser_update_wetdry_target(FFTDenoiser *self,
                                              const bool enable) {
  if (enable) {
    self->wet_dry_target = 1.f;
  } else {
    self->wet_dry_target = 0.f;
  }

  self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry) + FLT_MIN;
}

static void fft_denoiser_soft_bypass(FFTDenoiser *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->processed_fft_spectrum[k] =
        (1.f - self->wet_dry) * self->fft_spectrum[k] +
        self->processed_fft_spectrum[k] * self->wet_dry;
  }
}

static void residual_spectrum_whitening(FFTDenoiser *self,
                                        const float whitening_factor) {
  self->whitening_window_count++;

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->whitening_window_count > 1.f) {
      self->residual_max_spectrum[k] =
          fmaxf(fmaxf(self->residual_spectrum[k], WHITENING_FLOOR),
                self->residual_max_spectrum[k] * self->max_decay_rate);
    } else {
      self->residual_max_spectrum[k] =
          fmaxf(self->residual_spectrum[k], WHITENING_FLOOR);
    }
  }

  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    if (self->residual_spectrum[k] > FLT_MIN) {
      self->whitened_residual_spectrum[k] =
          self->residual_spectrum[k] / self->residual_max_spectrum[k];

      self->residual_spectrum[k] =
          (1.f - whitening_factor) * self->residual_spectrum[k] +
          whitening_factor * self->whitened_residual_spectrum[k];
    }
  }
}

static void get_denoised_spectrum(FFTDenoiser *self) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->denoised_spectrum[k] = self->fft_spectrum[k] * self->gain_spectrum[k];
  }
}

static void get_residual_spectrum(FFTDenoiser *self,
                                  const float whitening_factor) {
  for (uint32_t k = 1; k <= self->half_fft_size; k++) {
    self->residual_spectrum[k] =
        self->fft_spectrum[k] - self->denoised_spectrum[k];
  }

  if (whitening_factor > 0.f) {
    residual_spectrum_whitening(self, whitening_factor);
  }
}

static void get_final_spectrum(FFTDenoiser *self, const bool residual_listen,
                               const float reduction_amount) {
  if (residual_listen) {
    for (uint32_t k = 1; k <= self->half_fft_size; k++) {
      self->processed_fft_spectrum[k] = self->residual_spectrum[k];
    }
  } else {
    for (uint32_t k = 1; k <= self->half_fft_size; k++) {
      self->processed_fft_spectrum[k] =
          self->denoised_spectrum[k] +
          self->residual_spectrum[k] * reduction_amount;
    }
  }
}

static inline float from_db_to_coefficient(const float gain_db) {
  return expf(gain_db / 10.f * logf(10.f));
}
