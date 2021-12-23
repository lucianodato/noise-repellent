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

#include "noise_profile.h"
#include "string.h"
#include <stdlib.h>

struct NoiseProfile {
  uint32_t noise_profile_size;
  float *noise_profile;
};

NoiseProfile *noise_profile_initialize(const uint32_t size) {
  NoiseProfile *self = (NoiseProfile *)calloc(1U, sizeof(NoiseProfile));
  self->noise_profile_size = size;
  self->noise_profile = (float *)calloc(size, sizeof(float));

  return self;
}
void noise_profile_free(NoiseProfile *self) {
  free(self->noise_profile);
  free(self);
}

float *get_noise_profile(NoiseProfile *self) { return self->noise_profile; }

uint32_t get_noise_profile_size(NoiseProfile *self) {
  return self->noise_profile_size;
}

bool set_noise_profile(NoiseProfile *self, const float *noise_profile,
                       const uint32_t noise_profile_size) {
  if (!self || !noise_profile ||
      noise_profile_size != self->noise_profile_size) {
    return false;
  }
  memcpy(self->noise_profile, noise_profile,
         noise_profile_size * sizeof(float));

  return true;
}