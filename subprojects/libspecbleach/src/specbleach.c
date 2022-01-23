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

#include "../include/specbleach.h"
#include "adaptivedenoiser/adaptive_denoiser.h"
#include "denoiser/spectral_denoiser.h"
#include "shared/configurations.h"
#include "shared/general_utils.h"
#include "shared/noise_profile.h"
#include "stft/stft_processor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct SbSpectralDenoiser {
  uint32_t sample_rate;
  DenoiserParameters denoise_parameters;

  NoiseProfile *noise_profile;
  SpectralProcessorHandle spectral_denoiser;
  StftProcessor *stft_processor;
} SbSpectralDenoiser;

SpectralBleachHandle specbleach_initialize(const uint32_t sample_rate) {
  SbSpectralDenoiser *self =
      (SbSpectralDenoiser *)calloc(1U, sizeof(SbSpectralDenoiser));

  self->sample_rate = sample_rate;

  self->stft_processor = stft_processor_initialize(
      sample_rate, FRAME_SIZE_GENERAL, OVERLAP_FACTOR_GENERAL,
      PADDING_CONFIGURATION_GENERAL, INPUT_WINDOW_TYPE_GENERAL,
      OUTPUT_WINDOW_TYPE_GENERAL);

  if (!self->stft_processor) {
    specbleach_free(self);
    return NULL;
  }

  const uint32_t fft_size = get_stft_fft_size(self->stft_processor);
  const uint32_t real_spectrum_size =
      get_stft_real_spectrum_size(self->stft_processor);

  self->noise_profile = noise_profile_initialize(real_spectrum_size);

  if (!self->noise_profile) {
    specbleach_free(self);
    return NULL;
  }

  self->spectral_denoiser = spectral_denoiser_initialize(
      self->sample_rate, fft_size, OVERLAP_FACTOR_GENERAL, self->noise_profile);

  if (!self->spectral_denoiser) {
    specbleach_free(self);
    return NULL;
  }

  return self;
}

void specbleach_free(SpectralBleachHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  noise_profile_free(self->noise_profile);
  spectral_denoiser_free(self->spectral_denoiser);
  stft_processor_free(self->stft_processor);
  free(self);
}

uint32_t specbleach_get_latency(SpectralBleachHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  return get_stft_latency(self->stft_processor);
}

bool specbleach_process(SpectralBleachHandle instance,
                        const uint32_t number_of_samples, const float *input,
                        float *output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  stft_processor_run(self->stft_processor, number_of_samples, input, output,
                     &spectral_denoiser_run, self->spectral_denoiser);

  return true;
}

uint32_t specbleach_get_noise_profile_size(SpectralBleachHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  return get_noise_profile_size(self->noise_profile);
}

uint32_t
specbleach_get_noise_profile_blocks_averaged(SpectralBleachHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  return get_noise_profile_blocks_averaged(self->noise_profile);
}

float *specbleach_get_noise_profile(SpectralBleachHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  return get_noise_profile(self->noise_profile);
}

bool specbleach_load_noise_profile(SpectralBleachHandle instance,
                                   const float *restored_profile,
                                   const uint32_t profile_size,
                                   const uint32_t averaged_blocks) {
  if (!instance || !restored_profile) {
    return false;
  }

  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  if (profile_size != get_noise_profile_size(self->noise_profile)) {
    return false;
  }

  set_noise_profile(self->noise_profile, restored_profile, profile_size,
                    averaged_blocks);

  return true;
}

bool specbleach_reset_noise_profile(SpectralBleachHandle instance) {
  if (!instance) {
    return false;
  }

  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  reset_noise_profile(self->noise_profile);

  return true;
}

bool specbleach_noise_profile_available(SpectralBleachHandle instance) {
  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

  return is_noise_estimation_available(self->noise_profile);
}

bool specbleach_load_parameters(SpectralBleachHandle instance,
                                SpectralBleachParameters parameters) {
  if (!instance) {
    return false;
  }

  SbSpectralDenoiser *self = (SbSpectralDenoiser *)instance;

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

typedef struct SbAdaptiveDenoiser {
  uint32_t sample_rate;
  AdaptiveDenoiserParameters denoise_parameters;

  SpectralProcessorHandle adaptive_spectral_denoiser;
  StftProcessor *stft_processor;
} SbAdaptiveDenoiser;

SpectralBleachHandle
specbleach_adaptive_initialize(const uint32_t sample_rate) {
  SbAdaptiveDenoiser *self =
      (SbAdaptiveDenoiser *)calloc(1U, sizeof(SbAdaptiveDenoiser));

  self->sample_rate = sample_rate;

  self->stft_processor = stft_processor_initialize(
      sample_rate, FRAME_SIZE_SPEECH, OVERLAP_FACTOR_SPEECH,
      PADDING_CONFIGURATION_SPEECH, INPUT_WINDOW_TYPE_SPEECH,
      OUTPUT_WINDOW_TYPE_SPEECH);

  if (!self->stft_processor) {
    specbleach_free(self);
    return NULL;
  }

  const uint32_t fft_size = get_stft_fft_size(self->stft_processor);

  self->adaptive_spectral_denoiser =
      spectral_adaptive_denoiser_initialize(self->sample_rate, fft_size);

  if (!self->adaptive_spectral_denoiser) {
    specbleach_free(self);
    return NULL;
  }

  return self;
}

void specbleach_adaptive_free(SpectralBleachHandle instance) {
  SbAdaptiveDenoiser *self = (SbAdaptiveDenoiser *)instance;

  spectral_adaptive_denoiser_free(self->adaptive_spectral_denoiser);
  stft_processor_free(self->stft_processor);
  free(self);
}

uint32_t specbleach_adaptive_get_latency(SpectralBleachHandle instance) {
  SbAdaptiveDenoiser *self = (SbAdaptiveDenoiser *)instance;

  return get_stft_latency(self->stft_processor);
}

bool specbleach_adaptive_process(SpectralBleachHandle instance,
                                 const uint32_t number_of_samples,
                                 const float *input, float *output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  SbAdaptiveDenoiser *self = (SbAdaptiveDenoiser *)instance;

  stft_processor_run(self->stft_processor, number_of_samples, input, output,
                     &spectral_adaptive_denoiser_run,
                     self->adaptive_spectral_denoiser);

  return true;
}

bool specbleach_adaptive_load_parameters(SpectralBleachHandle instance,
                                         SpectralBleachParameters parameters) {
  if (!instance) {
    return false;
  }

  SbAdaptiveDenoiser *self = (SbAdaptiveDenoiser *)instance;

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