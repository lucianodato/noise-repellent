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

/**
* \file noise_profile.h
* \author Luciano Dato
* \brief The plugin state abstraction
*/

#ifndef NOISE_PROFILE_H
#define NOISE_PROFILE_H

typedef struct NoiseProfile NoiseProfile;

void set_noise_profile(NoiseProfile *self, float *noise_profile);
float *get_noise_profile(NoiseProfile *self);
void noise_profile_reset(NoiseProfile *self);
NoiseProfile *noise_profile_initialize(int noise_profile_size);
void noise_profile_free(NoiseProfile *self);

#endif