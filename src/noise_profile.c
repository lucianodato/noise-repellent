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
* \file noise_profile.c
* \author Luciano Dato
* \brief The plugin state abstraction
*/

#include "noise_profile.h"

struct NoiseProfile
{
	uint32_t child_size;
	uint32_t child_type;
	int noise_profile_size;
	float *values;
};

NoiseProfile *
noise_profile_initialize(LV2_URID child_type, int noise_profile_size)
{
	//Allocate object
	NoiseProfile *self = (NoiseProfile *)malloc(sizeof(NoiseProfile));

	self->child_type = child_type;
	self->child_size = sizeof(float);
	self->noise_profile_size = noise_profile_size;
	self->values = (float *)calloc((self->noise_profile_size), sizeof(float));

	return self;
}
