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

#ifndef SIGNAL_CROSSFADE_H
#define SIGNAL_CROSSFADE_H

#include <stdbool.h>
#include <stdint.h>

typedef struct SignalCrossfade SignalCrossfade;

SignalCrossfade* signal_crossfade_initialize(uint32_t sample_rate,
                                             uint32_t latency);
void signal_crossfade_free(SignalCrossfade* self);
bool signal_crossfade_run(SignalCrossfade* self, uint32_t number_of_samples,
                          const float* input, float* output, bool enable);
#endif
