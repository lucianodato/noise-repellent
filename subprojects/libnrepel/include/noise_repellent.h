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

// TODO -> Document interface

typedef struct NoiseRepellent NoiseRepellent;

NoiseRepellent *nr_initialize(uint32_t sample_rate);
void nr_free(NoiseRepellent *self);
bool nr_process(NoiseRepellent *self, uint32_t number_of_samples,
                const float *input, float *output);
uint32_t nr_get_latency(NoiseRepellent *self);
uint32_t nr_get_noise_profile_size(NoiseRepellent *self);
float *nr_get_noise_profile(NoiseRepellent *self);
bool nr_load_noise_profile(NoiseRepellent *self, const float *restored_profile);
bool nr_load_parameters(NoiseRepellent *self, bool enable, bool learn_noise,
                        float masking_ceiling_limit, float noise_rescale,
                        float reduction_amount, float release_time,
                        float residual_listen, float transient_threshold,
                        float whitening_factor); // TODO Expose parameters in
                                                 // struct and internal options

#endif