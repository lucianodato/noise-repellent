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

#include "../include/nrepel.h"
#include "adaptivedenoiser/adaptive_denoiser.h"
#include "denoiser/spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/general_utils.h"
#include "shared/noise_profile.h"
#include "stft/stft_processor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct NoiseRepellent {
  uint32_t sample_rate;
  DenoiserParameters denoise_parameters;

  NoiseProfile *noise_profile;
  SpectralProcessorHandle spectral_denoiser;
  StftProcessor *stft_processor;
} NoiseRepellent;

NoiseRepellentHandle nrepel_initialize(const uint32_t sample_rate) {
  NoiseRepellent *self = (NoiseRepellent *)calloc(1U, sizeof(NoiseRepellent));

  self->sample_rate = sample_rate;

  self->stft_processor = stft_processor_initialize(
      sample_rate, FRAME_SIZE_GENERAL, OVERLAP_FACTOR_GENERAL,
      INPUT_WINDOW_TYPE_GENERAL, OUTPUT_WINDOW_TYPE_GENERAL);

  if (!self->stft_processor) {
    nrepel_free(self);
    return NULL;
  }

  const uint32_t buffer_size = get_buffer_size(self->stft_processor);
  const uint32_t spectral_size =
      get_spectral_processing_size(self->stft_processor);

  self->noise_profile = noise_profile_initialize(spectral_size);

  if (!self->noise_profile) {
    nrepel_free(self);
    return NULL;
  }

  self->spectral_denoiser =
      spectral_denoiser_initialize(self->sample_rate, buffer_size,
                                   OVERLAP_FACTOR_GENERAL, self->noise_profile);

  if (!self->spectral_denoiser) {
    nrepel_free(self);
    return NULL;
  }

  return self;
}

void nrepel_free(NoiseRepellentHandle instance) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  noise_profile_free(self->noise_profile);
  spectral_denoiser_free(self->spectral_denoiser);
  stft_processor_free(self->stft_processor);
  free(self);
}

uint32_t nrepel_get_latency(NoiseRepellentHandle instance) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  return get_stft_latency(self->stft_processor);
}

bool nrepel_process(NoiseRepellentHandle instance,
                    const uint32_t number_of_samples, const float *input,
                    float *output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  NoiseRepellent *self = (NoiseRepellent *)instance;

  stft_processor_run(self->stft_processor, number_of_samples, input, output,
                     &spectral_denoiser_run, self->spectral_denoiser);

  return true;
}

uint32_t nrepel_get_noise_profile_size(NoiseRepellentHandle instance) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  return get_noise_profile_size(self->noise_profile);
}

uint32_t
nrepel_get_noise_profile_blocks_averaged(NoiseRepellentHandle instance) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  return get_noise_profile_blocks_averaged(self->noise_profile);
}

float *nrepel_get_noise_profile(NoiseRepellentHandle instance) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  return get_noise_profile(self->noise_profile);
}

bool nrepel_load_noise_profile(NoiseRepellentHandle instance,
                               const float *restored_profile,
                               const uint32_t profile_size,
                               const uint32_t averaged_blocks) {
  if (!instance || !restored_profile) {
    return false;
  }

  NoiseRepellent *self = (NoiseRepellent *)instance;

  if (profile_size != get_noise_profile_size(self->noise_profile)) {
    return false;
  }

  set_noise_profile(self->noise_profile, restored_profile, profile_size,
                    averaged_blocks);

  return true;
}

bool nrepel_reset_noise_profile(NoiseRepellentHandle instance) {
  if (!instance) {
    return false;
  }

  NoiseRepellent *self = (NoiseRepellent *)instance;

  reset_noise_profile(self->noise_profile);

  return true;
}

bool nrepel_noise_profile_available(NoiseRepellentHandle instance) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  return is_noise_estimation_available(self->noise_profile);
}

bool nrepel_load_parameters(NoiseRepellentHandle instance,
                            NrepelDenoiseParameters parameters) {
  if (!instance) {
    return false;
  }

  NoiseRepellent *self = (NoiseRepellent *)instance;

  // clang-format off
  self->denoise_parameters = (DenoiserParameters){
      .learn_noise = parameters.learn_noise,
      .residual_listen = parameters.residual_listen,
      .masking_ceiling_limit = parameters.masking_ceiling_limit,
      .reduction_amount =
          from_db_to_coefficient(parameters.reduction_amount * -1.F),
      .noise_rescale = from_db_to_coefficient(parameters.noise_rescale),
      .release_time = parameters.reduction_amount,
      .transient_threshold = parameters.transient_threshold,
      .whitening_factor = parameters.whitening_factor / 100.F,
  };
  // clang-format on

  load_reduction_parameters(self->spectral_denoiser, self->denoise_parameters);

  return true;
}

typedef struct NoiseRepellentAdaptive {
  uint32_t sample_rate;
  AdaptiveDenoiserParameters denoise_parameters;

  SpectralProcessorHandle adaptive_spectral_denoiser;
  StftProcessor *stft_processor;
} NoiseRepellentAdaptive;

NoiseRepellentHandle nrepel_adaptive_initialize(const uint32_t sample_rate) {
  NoiseRepellentAdaptive *self =
      (NoiseRepellentAdaptive *)calloc(1U, sizeof(NoiseRepellentAdaptive));

  self->sample_rate = sample_rate;

  self->stft_processor = stft_processor_initialize(
      sample_rate, FRAME_SIZE_SPEECH, OVERLAP_FACTOR_SPEECH,
      INPUT_WINDOW_TYPE_SPEECH, OUTPUT_WINDOW_TYPE_SPEECH);

  if (!self->stft_processor) {
    nrepel_free(self);
    return NULL;
  }

  const uint32_t buffer_size = get_buffer_size(self->stft_processor);

  self->adaptive_spectral_denoiser =
      spectral_adaptive_denoiser_initialize(self->sample_rate, buffer_size);

  if (!self->adaptive_spectral_denoiser) {
    nrepel_free(self);
    return NULL;
  }

  return self;
}

void nrepel_adaptive_free(NoiseRepellentHandle instance) {
  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  spectral_adaptive_denoiser_free(self->adaptive_spectral_denoiser);
  stft_processor_free(self->stft_processor);
  free(self);
}

uint32_t nrepel_adaptive_get_latency(NoiseRepellentHandle instance) {
  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  return get_stft_latency(self->stft_processor);
}

bool nrepel_adaptive_process(NoiseRepellentHandle instance,
                             const uint32_t number_of_samples,
                             const float *input, float *output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  stft_processor_run(self->stft_processor, number_of_samples, input, output,
                     &spectral_adaptive_denoiser_run,
                     self->adaptive_spectral_denoiser);

  return true;
}

bool nrepel_adaptive_load_parameters(NoiseRepellentHandle instance,
                                     NrepelDenoiseParameters parameters) {
  if (!instance) {
    return false;
  }

  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  // clang-format off
  self->denoise_parameters = (AdaptiveDenoiserParameters){
      .residual_listen = parameters.residual_listen,
      .reduction_amount =
          from_db_to_coefficient(parameters.reduction_amount * -1.F),
      .noise_rescale = from_db_to_coefficient(parameters.noise_rescale),
  };
  // clang-format on

  load_adaptive_reduction_parameters(self->adaptive_spectral_denoiser,
                                     self->denoise_parameters);

  return true;
}