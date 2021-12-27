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

#ifndef NREPEL_H
#define NREPEL_H

#include <stdbool.h>
#include <stdint.h>

// TODO (luciano/todo): Document interface
// TODO (luciano/todo): Test main file and increase coverage

typedef void *NoiseRepellentHandle;

typedef struct {
  bool enable;
  bool learn_noise;
  bool residual_listen;
  bool adaptive_noise_learn;
  float reduction_amount;
  float release_time;
  float masking_ceiling_limit;
  float whitening_factor;
  float transient_threshold;
  float noise_rescale;
} ProcessorParameters;

NoiseRepellentHandle nr_initialize(uint32_t sample_rate);
void nr_free(NoiseRepellentHandle instance);
bool nr_process(NoiseRepellentHandle instance, uint32_t number_of_samples,
                const float *input, float *output);
uint32_t nr_get_latency(NoiseRepellentHandle instance);
uint32_t nr_get_noise_profile_size(NoiseRepellentHandle instance);
float *nr_get_noise_profile(NoiseRepellentHandle instance);
bool nr_load_noise_profile(NoiseRepellentHandle instance,
                           const float *restored_profile, uint32_t profile_size,
                           uint32_t profile_blocks);
uint32_t nr_get_noise_profile_blocks_averaged(NoiseRepellentHandle instance);
bool nr_reset_noise_profile(NoiseRepellentHandle instance);
bool nr_load_parameters(NoiseRepellentHandle instance,
                        ProcessorParameters parameters);
#endif