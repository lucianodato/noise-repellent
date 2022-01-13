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

#ifndef STFT_WINDOW_H
#define STFT_WINDOW_H

#include "../shared/spectral_utils.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct StftWindows StftWindows;

typedef enum WindowPlace { INPUT_WINDOW = 1, OUTPUT_WINDOW = 2 } WindowPlace;

StftWindows *stft_window_initialize(uint32_t window_size,
                                    uint32_t overlap_factor,
                                    WindowTypes input_window,
                                    WindowTypes output_window);
void stft_window_free(StftWindows *self);
bool apply_window(StftWindows *self, float *frame, WindowPlace place);

#endif