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

#include "../include/noise_repellent.h"
#include "denoiser/spectral_denoiser.h"
#include "noisemodel/noise_estimator.h"
#include "shared/common.h"
#include "shared/noise_profile.h"
#include "shared/signal_crossfade.h"
#include "stft/stft_processor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseRepellent {
  ProcessorParameters *denoise_parameters;
  NoiseProfile *noise_profile;

  NoiseEstimator *noise_estimator;
  SpectralDenoiser *spectral_denoiser;
  StftProcessor *stft_processor;
  SignalCrossfade *soft_bypass;

  uint32_t sample_rate;
};

NoiseRepellent *nr_initialize(const uint32_t sample_rate) {
  NoiseRepellent *self = (NoiseRepellent *)calloc(1, sizeof(NoiseRepellent));
  self->sample_rate = sample_rate;

  self->stft_processor = stft_processor_initialize();

  if (!self->stft_processor) {
    nr_free(self);
    return NULL;
  }

  const uint32_t fft_size = get_fft_size(self->stft_processor);
  const uint32_t overlap_factor = get_overlap_factor(self->stft_processor);
  const uint32_t spectral_size =
      get_spectral_processing_size(self->stft_processor);

  self->noise_profile = noise_profile_initialize(spectral_size);

  if (!self->noise_profile) {
    nr_free(self);
    return NULL;
  }

  self->soft_bypass = signal_crossfade_initialize(self->sample_rate);

  if (!self->soft_bypass) {
    nr_free(self);
    return NULL;
  }

  self->denoise_parameters =
      (ProcessorParameters *)calloc(1, sizeof(ProcessorParameters));

  if (!self->denoise_parameters) {
    nr_free(self);
    return NULL;
  }

  self->noise_estimator =
      noise_estimation_initialize(fft_size, self->noise_profile);

  if (!self->noise_estimator) {
    nr_free(self);
    return NULL;
  }

  self->spectral_denoiser = spectral_denoiser_initialize(
      self->sample_rate, fft_size, overlap_factor, self->noise_profile,
      self->denoise_parameters);

  if (!self->spectral_denoiser) {
    nr_free(self);
    return NULL;
  }

  return self;
}

void nr_free(NoiseRepellent *self) {
  free(self->denoise_parameters);
  signal_crossfade_free(self->soft_bypass);
  noise_profile_free(self->noise_profile);
  noise_estimation_free(self->noise_estimator);
  spectral_denoiser_free(self->spectral_denoiser);
  stft_processor_free(self->stft_processor);
  free(self);
}

uint32_t nr_get_latency(NoiseRepellent *self) {
  return get_stft_latency(self->stft_processor);
}

bool nr_process(NoiseRepellent *self, const uint32_t number_of_samples,
                const float *input, float *output) {
  if (!self || !number_of_samples || !input || !output) {
    return false;
  }

  if (self->denoise_parameters->learn_noise) {
    stft_processor_run(self->stft_processor, &noise_estimation_run,
                       (SPECTRAL_PROCESSOR)self->noise_estimator,
                       number_of_samples, input, output); // estimating noise
  } else if (is_noise_estimation_available(self->noise_estimator)) {
    stft_processor_run(self->stft_processor, &spectral_denoiser_run,
                       (SPECTRAL_PROCESSOR)self->spectral_denoiser,
                       number_of_samples, input, output); // denoising
  } else {
    memcpy(output, input,
           sizeof(float) * number_of_samples); // bypassed
  }

  signal_crossfade_run(self->soft_bypass, number_of_samples, input, output,
                       self->denoise_parameters->enable);

  return true;
}

uint32_t nr_get_noise_profile_size(NoiseRepellent *self) {
  return get_noise_profile_size(self->noise_profile);
}

float *nr_get_noise_profile(NoiseRepellent *self) {
  return get_noise_profile(self->noise_profile);
}

bool nr_load_noise_profile(NoiseRepellent *self,
                           const float *restored_profile) {
  if (!self || !restored_profile) { // TODO check if size corresponds to the
                                    // configured one
    return false;
  }

  set_noise_profile(self->noise_profile, restored_profile,
                    get_noise_profile_size(self->noise_profile));

  return true;
}

static inline float from_db_to_coefficient(const float gain_db) {
  return expf(gain_db / 10.f * logf(10.f));
}

bool nr_load_parameters(NoiseRepellent *self, const bool enable,
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

  return true;
}