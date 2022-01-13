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
#include "shared/configurations.h"
#include "shared/general_utils.h"
#include "stft/stft_processor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct NoiseRepellentAdaptive {
  uint32_t sample_rate;
  NrepelDenoiseParameters denoise_parameters;

  SpectralProcessorHandle adaptive_spectral_denoiser;
  StftProcessor *stft_processor;
} NoiseRepellentAdaptive;

NoiseRepellentHandle nrepel_initialize(const uint32_t sample_rate) {
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

void nrepel_free(NoiseRepellentHandle instance) {
  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  spectral_adaptive_denoiser_free(self->adaptive_spectral_denoiser);
  stft_processor_free(self->stft_processor);
  free(self);
}

uint32_t nrepel_get_latency(NoiseRepellentHandle instance) {
  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  return get_stft_latency(self->stft_processor);
}

bool nrepel_process(NoiseRepellentHandle instance,
                    const uint32_t number_of_samples, const float *input,
                    float *output) {
  if (!instance || number_of_samples == 0 || !input || !output) {
    return false;
  }

  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  load_adaptive_reduction_parameters(self->adaptive_spectral_denoiser,
                                     self->denoise_parameters.residual_listen,
                                     self->denoise_parameters.reduction_amount,
                                     self->denoise_parameters.noise_rescale);
  stft_processor_run(self->stft_processor, number_of_samples, input, output,
                     &spectral_adaptive_denoiser_run,
                     self->adaptive_spectral_denoiser);

  return true;
}

bool nrepel_load_parameters(NoiseRepellentHandle instance,
                            NrepelDenoiseParameters parameters) {
  if (!instance) {
    return false;
  }

  NoiseRepellentAdaptive *self = (NoiseRepellentAdaptive *)instance;

  // clang-format off
  self->denoise_parameters = (NrepelDenoiseParameters){
      .residual_listen = parameters.residual_listen,
      .reduction_amount =
          from_db_to_coefficient(parameters.reduction_amount * -1.F),
      .noise_rescale = from_db_to_coefficient(parameters.noise_rescale),
  };
  // clang-format on

  return true;
}