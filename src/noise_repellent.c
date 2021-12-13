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

#include "noise_repellent.h"
#include "fft_denoiser.h"
#include "stft_processor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NOISE_PROFILE_SIZE FFT_SIZE / 2 + 1

struct NoiseRepellentLib {
  DenoiseParameters *denoise_parameters;
  NoiseProfile *noise_profile;

  FFTDenoiser *fft_denoiser;
  STFTProcessor *stft_processor;

  uint32_t sample_rate;
};

NoiseRepellentLib *nr_initialize(const uint32_t sample_rate) {
  NoiseRepellentLib *self =
      (NoiseRepellentLib *)calloc(1, sizeof(NoiseRepellentLib));
  self->sample_rate = sample_rate;

  self->denoise_parameters =
      (DenoiseParameters *)calloc(1, sizeof(DenoiseParameters));

  self->noise_profile = (NoiseProfile *)calloc(1, sizeof(NoiseProfile));
  self->noise_profile->noise_profile_size = NOISE_PROFILE_SIZE;
  self->noise_profile->noise_profile =
      (float *)calloc(NOISE_PROFILE_SIZE, sizeof(float));

  self->fft_denoiser =
      fft_denoiser_initialize(self->sample_rate, FFT_SIZE, OVERLAP_FACTOR);
  self->stft_processor = stft_processor_initialize(FFT_SIZE, OVERLAP_FACTOR);

  if (!self->denoise_parameters || !self->noise_profile ||
      !self->noise_profile->noise_profile || !self->fft_denoiser ||
      !self->stft_processor) {
    nr_free(self);
    return NULL;
  }

  load_noise_profile(self->fft_denoiser, self->noise_profile);
  load_denoiser(self->stft_processor, self->fft_denoiser);

  return self;
}

void nr_free(NoiseRepellentLib *self) {
  free(self->noise_profile->noise_profile);
  free(self->noise_profile);
  free(self->denoise_parameters);
  fft_denoiser_free(self->fft_denoiser);
  stft_processor_free(self->stft_processor);
  free(self);
}

uint32_t nr_get_latency(NoiseRepellentLib *self) {
  return get_stft_latency(self->stft_processor);
}

bool nr_process(NoiseRepellentLib *self, const uint32_t number_of_samples,
                const float *input, float *output) {
  if (!self || !number_of_samples || !input || !output) {
    return false;
  }
  stft_processor_run(self->stft_processor, number_of_samples, input, output);
  return true;
}

uint32_t nr_get_noise_profile_size(NoiseRepellentLib *self) {
  return NOISE_PROFILE_SIZE;
}

float *nr_get_noise_profile(NoiseRepellentLib *self) {
  return self->noise_profile->noise_profile;
}

bool nr_load_noise_profile(NoiseRepellentLib *self,
                           const float restored_profile[]) {
  if (!self || !restored_profile) {
    return false;
  }
  memcpy(self->noise_profile->noise_profile, restored_profile,
         sizeof(float) * NOISE_PROFILE_SIZE);
  return true;
}

static inline float from_db_to_coefficient(const float gain_db) {
  return expf(gain_db / 10.f * logf(10.f));
}

bool nr_load_parameters(NoiseRepellentLib *self, const bool enable,
                        const bool learn_noise,
                        const float masking_ceiling_limit,
                        const float noise_rescale, const float reduction_amount,
                        const float release_time, const float residual_listen,
                        const float transient_threshold,
                        const float whitening_factor) {
  if (!self) {
    return false;
  }

  self->denoise_parameters->enable = enable;
  self->denoise_parameters->learn_noise = learn_noise;
  self->denoise_parameters->residual_listen = residual_listen;
  self->denoise_parameters->reduction_amount =
      from_db_to_coefficient(reduction_amount * -1.f);
  self->denoise_parameters->release_time = release_time;
  self->denoise_parameters->masking_ceiling_limit =
      masking_ceiling_limit / 100.f;
  self->denoise_parameters->whitening_factor = whitening_factor;
  self->denoise_parameters->transient_threshold = transient_threshold;
  self->denoise_parameters->noise_rescale = noise_rescale;

  load_denoise_parameters(self->fft_denoiser, self->denoise_parameters);

  return true;
}