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
* \file plugin_state.h
* \author Luciano Dato
* \brief The plugin state abstraction
*/

#ifndef PLUGIN_STATE_H
#define PLUGIN_STATE_H

#include "lv2/lv2plug.in/ns/ext/atom/atom.h"
#include "lv2/lv2plug.in/ns/ext/state/state.h"
#include "lv2/lv2plug.in/ns/ext/urid/urid.h"
#include "lv2/lv2plug.in/ns/lv2core/lv2.h"
#include "noise_profile.h"

#define NOISEREPELLENT_URI "https://github.com/lucianodato/noise-repellent"

typedef struct PluginState PluginState;

bool plugin_state_initialize(PluginState *self, const LV2_Feature *const *features);
void plugin_state_free(PluginState *self);
void plugin_state_savestate(PluginState *self, LV2_State_Store_Function store, LV2_State_Handle handle,
							int fft_size, NoiseProfile *noise_profile);
bool plugin_state_restorestate(PluginState *self, LV2_State_Retrieve_Function retrieve, LV2_State_Handle handle,
							   NoiseProfile *noise_profile, int *fft_size);

#endif