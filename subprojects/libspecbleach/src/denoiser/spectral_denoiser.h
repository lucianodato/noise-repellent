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

#ifndef SPECTRAL_DENOISER_H
#define SPECTRAL_DENOISER_H

#include "../shared/noise_profile.h"
#include "../shared/spectral_processor.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct DenoiserParameters {
  float reduction_amount;
  float noise_rescale;
  bool residual_listen;
  bool learn_noise;
  float release_time;
  float whitening_factor;
  float transient_threshold;
} DenoiserParameters;

SpectralProcessorHandle
spectral_denoiser_initialize(uint32_t sample_rate, uint32_t fft_size,
                             uint32_t overlap_factor,
                             NoiseProfile *noise_profile);
void spectral_denoiser_free(SpectralProcessorHandle instance);
bool load_reduction_parameters(SpectralProcessorHandle instance,
                               DenoiserParameters parameters);
bool spectral_denoiser_run(SpectralProcessorHandle instance,
                           float *fft_spectrum);

#endif