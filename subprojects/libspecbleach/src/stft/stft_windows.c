/*
libspecbleach - A spectral processing library

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

#include "stft_windows.h"
#include "../shared/configurations.h"
#include <stdlib.h>

static float get_windows_scale_factor(StftWindows *self,
                                      uint32_t overlap_factor);

struct StftWindows {
  float *input_window;
  float *output_window;
  uint32_t stft_frame_size;
  float scale_factor;
};

StftWindows *stft_window_initialize(const uint32_t stft_frame_size,
                                    const uint32_t overlap_factor,
                                    const WindowTypes input_window,
                                    const WindowTypes output_window) {
  StftWindows *self = (StftWindows *)calloc(1U, sizeof(StftWindows));

  self->stft_frame_size = stft_frame_size;

  self->input_window = (float *)calloc(self->stft_frame_size, sizeof(float));
  self->output_window = (float *)calloc(self->stft_frame_size, sizeof(float));

  get_fft_window(self->input_window, self->stft_frame_size, input_window);
  get_fft_window(self->output_window, self->stft_frame_size, output_window);

  self->scale_factor = get_windows_scale_factor(self, overlap_factor);

  return self;
}

void stft_window_free(StftWindows *self) {
  free(self->input_window);
  free(self->output_window);

  free(self);
}

static float get_windows_scale_factor(StftWindows *self,
                                      const uint32_t overlap_factor) {
  if (overlap_factor < 2) {
    return 0.F;
  }
  float sum = 0.F;
  for (uint32_t i = 0U; i < self->stft_frame_size; i++) {
    sum += self->input_window[i] * self->output_window[i];
  }

  return sum * (float)overlap_factor;
}

bool stft_window_apply(StftWindows *self, float *frame,
                       const WindowPlace place) {
  if (!self || !frame) {
    return false;
  }

  for (uint32_t i = 0U; i < self->stft_frame_size; i++) {
    switch (place) {
    case INPUT_WINDOW:
      frame[i] *= self->input_window[i];
      break;
    case OUTPUT_WINDOW:
      frame[i] *= self->output_window[i] / self->scale_factor;
      break;
    default:
      break;
    }
  }

  return true;
}