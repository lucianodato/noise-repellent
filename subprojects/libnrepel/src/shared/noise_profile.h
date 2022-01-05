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

#ifndef NOISE_PROFILE_H
#define NOISE_PROFILE_H

#include <stdbool.h>
#include <stdint.h>

typedef struct NoiseProfile NoiseProfile;

NoiseProfile *noise_profile_initialize(uint32_t size);
void noise_profile_free(NoiseProfile *self);
float *get_noise_profile(NoiseProfile *self);
uint32_t get_noise_profile_size(NoiseProfile *self);
uint32_t get_noise_profile_blocks_averaged(NoiseProfile *self);
bool increment_blocks_averaged(NoiseProfile *self);
bool set_noise_profile(NoiseProfile *self, const float *noise_profile,
                       uint32_t noise_profile_size, uint32_t averaged_blocks);
bool reset_noise_profile(NoiseProfile *self);
bool is_noise_estimation_available(NoiseProfile *self);
void set_noise_estimation_available(NoiseProfile *self);

#endif
