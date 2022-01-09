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

#include "stft_windows.h"
#include "../shared/configurations.h"
#include "../shared/spectral_utils.h"
#include <stdlib.h>

static float get_windows_scale_factor(StftWindows *self,
                                      uint32_t overlap_factor);

struct StftWindows {
  uint32_t window_option_input;
  uint32_t window_option_output;
  float *input_window;
  float *output_window;
  uint32_t window_size;
  float scale_factor;
};

StftWindows *stft_window_initialize(const uint32_t window_size,
                                    const uint32_t overlap_factor) {
  StftWindows *self = (StftWindows *)calloc(1U, sizeof(StftWindows));

  self->window_size = window_size;
  self->window_option_input = INPUT_WINDOW_TYPE;
  self->window_option_output = OUTPUT_WINDOW_TYPE;

  self->input_window = (float *)calloc(self->window_size, sizeof(float));
  self->output_window = (float *)calloc(self->window_size, sizeof(float));

  get_fft_window(self->input_window, self->window_size,
                 self->window_option_input);
  get_fft_window(self->output_window, self->window_size,
                 self->window_option_output);

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
  for (uint32_t i = 0U; i < self->window_size; i++) {
    sum += self->input_window[i] * self->output_window[i];
  }

  return sum * (float)overlap_factor;
}

bool apply_window(StftWindows *self, float *frame, const WindowPlace place) {
  if (!self || !frame) {
    return false;
  }

  for (uint32_t i = 0U; i < self->window_size / 2U; i++) {
    switch (place) {
    case INPUT_WINDOW:
      frame[i] *= self->input_window[i];
      frame[self->window_size - 1U - i] *=
          self->input_window[self->window_size - 1U - i];
      break;
    case OUTPUT_WINDOW:
      frame[i] *= self->output_window[i] / self->scale_factor;
      frame[self->window_size - 1U - i] *=
          self->output_window[self->window_size - 1U - i] / self->scale_factor;
      break;
    default:
      break;
    }
  }

  return true;
}