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

#ifndef NREPEL_H_INCLUDED
#define NREPEL_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>

// TODO (luciano/todo): Extract library to it's own repository when API is
// stable. Manage visibility with meson
// TODO (luciano/todo): Document interface
// TODO (luciano/todo): Test main file and increase coverage
// TODO (luciano/todo): Add post filter for frequency smoothing
// TODO (luciano/todo): Move Plugin functionality to the plugin.
// noise_profile instance, parameters, etc.

typedef void *NoiseRepellentHandle;

typedef struct NrepelDenoiseParameters {
  bool learn_noise;     // Plugin
  bool residual_listen; // Plugin
  float reduction_amount;
  float release_time;
  float masking_ceiling_limit; // Should be internal
  float whitening_factor;
  float transient_threshold; // Should be Adaptive or fixed
  float noise_rescale;
} NrepelDenoiseParameters;

NoiseRepellentHandle nrepel_initialize(uint32_t sample_rate);
void nrepel_free(NoiseRepellentHandle instance);
bool nrepel_process(NoiseRepellentHandle instance, uint32_t number_of_samples,
                    const float *input, float *output);
uint32_t nrepel_get_latency(NoiseRepellentHandle instance);
bool nrepel_load_parameters(NoiseRepellentHandle instance,
                            NrepelDenoiseParameters parameters);
uint32_t nrepel_get_noise_profile_size(NoiseRepellentHandle instance);
float *nrepel_get_noise_profile(NoiseRepellentHandle instance);
bool nrepel_load_noise_profile(NoiseRepellentHandle instance,
                               const float *restored_profile,
                               uint32_t profile_size, uint32_t profile_blocks);
uint32_t
nrepel_get_noise_profile_blocks_averaged(NoiseRepellentHandle instance);
bool nrepel_reset_noise_profile(NoiseRepellentHandle instance);
bool nrepel_noise_profile_available(NoiseRepellentHandle instance);

NoiseRepellentHandle nrepel_adaptive_initialize(uint32_t sample_rate);
void nrepel_adaptive_free(NoiseRepellentHandle instance);
bool nrepel_adaptive_process(NoiseRepellentHandle instance,
                             uint32_t number_of_samples, const float *input,
                             float *output);
uint32_t nrepel_adaptive_get_latency(NoiseRepellentHandle instance);
bool nrepel_adaptive_load_parameters(NoiseRepellentHandle instance,
                                     NrepelDenoiseParameters parameters);

#endif