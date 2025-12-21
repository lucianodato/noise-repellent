/*
noise-repellent -- Noise Reduction LV2

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "noise_profile_state.h"

#define MAX_PROFILE_SIZE 8192

struct NoiseProfileState {
  uint32_t child_size;
  uint32_t child_type;
  float elements[MAX_PROFILE_SIZE];
}; // LV2 Atoms Vector Specification

NoiseProfileState* noise_profile_state_initialize(LV2_URID child_type) {
  NoiseProfileState* self =
      (NoiseProfileState*)calloc(1U, sizeof(NoiseProfileState));
  self->child_type = (uint32_t)child_type;
  self->child_size = (uint32_t)sizeof(float);

  return self;
}

void noise_profile_state_free(NoiseProfileState* self) {
  free(self);
}

float* noise_profile_get_elements(NoiseProfileState* self) {
  return self->elements;
}
size_t noise_profile_get_size(void) {
  return sizeof(NoiseProfileState);
}
