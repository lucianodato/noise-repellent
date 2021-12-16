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
#include "data_types.h"
#include "spectral_processor.h"
#include "stft_processor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct NoiseRepellent {
  ProcessorParameters *denoise_parameters;
  NoiseProfile *noise_profile;

  SpectralProcessor *fft_denoiser;
  StftProcessor *stft_processor;

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

  self->fft_denoiser = spectral_processor_initialize(self->sample_rate,
                                                     fft_size, overlap_factor);

  if (!self->fft_denoiser) {
    nr_free(self);
    return NULL;
  }

  self->denoise_parameters =
      (ProcessorParameters *)calloc(1, sizeof(ProcessorParameters));

  if (!self->denoise_parameters) {
    nr_free(self);
    return NULL;
  }

  self->noise_profile = (NoiseProfile *)calloc(1, sizeof(NoiseProfile));
  self->noise_profile->noise_profile =
      (float *)calloc(spectral_size, sizeof(float));

  if (!self->noise_profile || !self->noise_profile->noise_profile) {
    nr_free(self);
    return NULL;
  }

  self->noise_profile->noise_profile_size = spectral_size;
  load_noise_profile(self->fft_denoiser, self->noise_profile);

  return self;
}

void nr_free(NoiseRepellent *self) {
  free(self->noise_profile);
  free(self->denoise_parameters);
  spectral_processor_free(self->fft_denoiser);
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

  stft_processor_run(self->stft_processor, &spectral_processor_run,
                     (SPECTAL_PROCESSOR)self->fft_denoiser, number_of_samples,
                     input, output);
  return true;
}

uint32_t nr_get_noise_profile_size(NoiseRepellent *self) {
  return self->noise_profile->noise_profile_size;
}

float *nr_get_noise_profile(NoiseRepellent *self) {
  return self->noise_profile->noise_profile;
}

bool nr_load_noise_profile(NoiseRepellent *self,
                           const float restored_profile[]) {
  if (!self || !restored_profile) {
    return false;
  }
  memcpy(self->noise_profile->noise_profile, restored_profile,
         sizeof(float) * self->noise_profile->noise_profile_size);
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

  load_processor_parameters(self->fft_denoiser, self->denoise_parameters);

  return true;
}