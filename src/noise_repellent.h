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

#ifndef NOISE_REPELLENT_H
#define NOISE_REPELLENT_H

#include <stdbool.h>
#include <stdint.h>

typedef struct NoiseRepellentLib NoiseRepellentLib;

NoiseRepellentLib *nr_initialize(uint32_t sample_rate);
void nr_free(NoiseRepellentLib *self);
uint32_t nr_get_latency(NoiseRepellentLib *self);
bool nr_process(NoiseRepellentLib *self, uint32_t number_of_samples,
                const float *input, float *output);
uint32_t nr_get_noise_profile_size(NoiseRepellentLib *self);
float *nr_get_noise_profile(NoiseRepellentLib *self);
bool nr_load_noise_profile(NoiseRepellentLib *self,
                           const float *restored_profile);
bool nr_load_parameters(NoiseRepellentLib *self, bool enable, bool learn_noise,
                        float masking_ceiling_limit, float noise_rescale,
                        float reduction_amount, float release_time,
                        float residual_listen, float transient_threshold,
                        float whitening_factor);

#endif