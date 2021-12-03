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
#include <stdlib.h>
#include <string.h>

struct NoiseProfile
{
	int noise_profile_size;
	float *values;
};

void set_noise_profile(NoiseProfile *self, float *noise_profile)
{
	if (*noise_profile)
	{
		memcpy(self->values, noise_profile, self->noise_profile_size);
	}
}

float *get_noise_profile(NoiseProfile *self)
{
	return self->values;
}

void noise_profile_reset(NoiseProfile *self)
{
	memset(self->values, 0.f, self->noise_profile_size + 1);
}

NoiseProfile *noise_profile_initialize(int noise_profile_size)
{
	//Allocate object
	NoiseProfile *self = (NoiseProfile *)malloc(sizeof(NoiseProfile));

	self->noise_profile_size = noise_profile_size;
	self->values = (float *)calloc((self->noise_profile_size), sizeof(float));

	noise_profile_reset(self);

	return self;
}

void noise_profile_free(NoiseProfile *self)
{
	free(self->values);
	free(self);
}
