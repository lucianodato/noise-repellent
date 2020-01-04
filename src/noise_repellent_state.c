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
* \file plugin_state.c
* \author Luciano Dato
* \brief The plugin state abstraction
*/

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"
#include "lv2/lv2plug.in/ns/ext/urid/urid.h"
#include "lv2/lv2plug.in/ns/ext/atom/atom.h"
#include "lv2/lv2plug.in/ns/ext/state/state.h"
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

#define NOISEREPELLENT_URI "https://github.com/lucianodato/noise-repellent"

/**
* Noise Profile state.
*/
typedef struct
{
    uint32_t child_size;
    uint32_t child_type;
    int np_size;
    float *values;
} NoiseProfile;

NoiseProfile*
np_init(LV2_URID child_type, int np_size)
{
    //Allocate object
    NoiseProfile *self = (NoiseProfile *)malloc(sizeof(NoiseProfile));

    self->child_type = child_type;
	self->child_size = sizeof(float);
    self->np_size = np_size;
    self->values = (float *)calloc((self->np_size), sizeof(float));

    return self;
}

/**
* Struct for the plugin state.
*/
typedef struct
{
    //LV2 state URID (Save and restore noise profile)
    LV2_URID_Map *map;
    LV2_URID atom_Vector;
    LV2_URID atom_Int;
    LV2_URID atom_Float;
    LV2_URID prop_fftsize;
    LV2_URID prop_nwindow;
    LV2_URID prop_FFTp2;

    NoiseProfile *noise_profile;
} Plugin_State;

bool ps_configure(Plugin_State *self, const LV2_Feature *const *features, int np_size)
{
    //Retrieve the URID map callback, and needed URIDs
    for (int i = 0; features[i]; ++i)
    {
        if (!strcmp(features[i]->URI, LV2_URID__map))
        {
            self->map = (LV2_URID_Map *)features[i]->data;
        }
    }
    if (!self->map)
    {
        return false; //host doesn't support states
    }

    //For lv2 state (noise profile saving)
    self->atom_Vector = self->map->map(self->map->handle, LV2_ATOM__Vector);
    self->atom_Int = self->map->map(self->map->handle, LV2_ATOM__Int);
    self->atom_Float = self->map->map(self->map->handle, LV2_ATOM__Float);
    self->prop_fftsize = self->map->map(self->map->handle, NOISEREPELLENT_URI "#fftsize");
    self->prop_nwindow = self->map->map(self->map->handle, NOISEREPELLENT_URI "#nwindow");
    self->prop_FFTp2 = self->map->map(self->map->handle, NOISEREPELLENT_URI "#FFTp2");

    self->noise_profile = np_init(self->atom_Float,np_size);

    return true;
}

void ps_savestate(Plugin_State *self, LV2_State_Store_Function store, LV2_State_Handle handle,
                  int *fft_size,float *noise_window_count,float *noise_profile)
{
    store(handle, self->prop_fftsize, &fft_size, sizeof(int), self->atom_Int,
        LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	store(handle, self->prop_nwindow, &noise_window_count, sizeof(float),
		  self->atom_Float, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	memcpy(self->noise_profile->values, noise_profile, sizeof(self->noise_profile->np_size));

	store(handle, self->prop_FFTp2, (void *)self->noise_profile, sizeof(NoiseProfile),
		  self->atom_Vector, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);
}

bool ps_restorestate(Plugin_State *self,  LV2_State_Retrieve_Function retrieve, LV2_State_Handle handle,
                     float* noise_profile, float noise_window_count, int *fft_size,
                     int fft_size_2)
{
    size_t size;
	uint32_t type;
	uint32_t valflags;

	const int32_t *fftsize = retrieve(handle, self->prop_fftsize, &size, &type, &valflags);
	if (!fftsize || type != self->atom_Int || *fftsize != fft_size_2)
	{
		return false;
	}

	const void *vecFFTp2 = retrieve(handle, self->prop_FFTp2, &size, &type, &valflags);
	if (!vecFFTp2 || size != sizeof(NoiseProfile) || type != self->atom_Vector)
	{
		return false;
	}

	//Deactivate any denoising before loading any noise profile
	//self->noise_thresholds_availables = false;

	//Copy to local variables
	memcpy(self->noise_profile, (float *)LV2_ATOM_BODY(vecFFTp2), (fft_size_2 + 1) * sizeof(float));

	const float *wincount = retrieve(handle, self->prop_nwindow, &size, &type, &valflags);
	if (fftsize && type == self->atom_Float)
	{
		noise_window_count = *wincount;
	}

	//Reactivate denoising with restored profile
	//self->noise_thresholds_availables = true;

    return true;
}