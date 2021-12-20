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

#include "spectral_utils.h"
#include "common.h"
#include <math.h>

static inline float blackman(const uint32_t bin_index,
                             const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return 0.42 - 0.5 * cosf(2.f * M_PI * p) + 0.08 * cosf(4.f * M_PI * p);
}

static inline float hanning(const uint32_t bin_index, const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return 0.5 - 0.5 * cosf(2.f * M_PI * p);
}

static inline float hamming(const uint32_t bin_index, const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return 0.54 - 0.46 * cosf(2.f * M_PI * p);
}

static inline float vorbis(const uint32_t bin_index, const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return sinf(M_PI / 2.f * powf(sinf(M_PI * p), 2.f));
}

bool get_fft_window(float *window, const uint32_t fft_size,
                    const WindowTypes window_type) {
  if (!window || !fft_size || !window_type) {
    return false;
  }

  for (uint32_t k = 0; k < fft_size; k++) {
    switch (window_type) {
    case BLACKMAN_WINDOW:
      window[k] = blackman(k, fft_size);
      break;
    case HANN_WINDOW:
      window[k] = hanning(k, fft_size);
      break;
    case HAMMING_WINDOW:
      window[k] = hamming(k, fft_size);
      break;
    case VORBIS_WINDOW:
      window[k] = vorbis(k, fft_size);
      break;
    }
  }

  return true;
}