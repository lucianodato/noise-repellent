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

#include "spectral_utils.h"
#include "configurations.h"
#include "general_utils.h"
#include <float.h>
#include <math.h>

static inline float blackman(const uint32_t bin_index,
                             const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return sanitize_denormal(0.42F - (0.5F * cosf(2.F * M_PI * p)) +
                           (0.08F * cosf(4.F * M_PI * p)));
}

static inline float hanning(const uint32_t bin_index, const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return sanitize_denormal(0.5F - (0.5F * cosf(2.F * M_PI * p)));
}

static inline float hamming(const uint32_t bin_index, const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return sanitize_denormal(0.54F - (0.46F * cosf(2.F * M_PI * p)));
}

static inline float vorbis(const uint32_t bin_index, const uint32_t fft_size) {
  const float p = ((float)(bin_index)) / ((float)(fft_size));
  return sanitize_denormal(sinf(M_PI / 2.F * powf(sinf(M_PI * p), 2.F)));
}

bool get_fft_window(float *window, const uint32_t fft_size,
                    const WindowTypes window_type) {
  if (!window || !fft_size) {
    return false;
  }

  for (uint32_t k = 0; k < fft_size; k++) {
    switch (window_type) {
    case HANN_WINDOW:
      window[k] = hanning(k, fft_size);
      break;
    case HAMMING_WINDOW:
      window[k] = hamming(k, fft_size);
      break;
    case BLACKMAN_WINDOW:
      window[k] = blackman(k, fft_size);
      break;
    case VORBIS_WINDOW:
      window[k] = vorbis(k, fft_size);
      break;
    default:
      break;
    }
  }

  return true;
}

bool initialize_spectrum_with_value(float *spectrum, uint32_t spectrum_size,
                                    const float value) {
  if (!spectrum || spectrum_size <= 0U) {
    return false;
  }

  for (uint32_t i = 0U; i < spectrum_size; i++) {
    spectrum[i] = value;
  }

  return true;
}

float max_spectral_value(const float *spectrum, const uint32_t real_spectrum_size) {
  if (!spectrum || real_spectrum_size <= 0U) {
    return 0.F;
  }

  float max = spectrum[0];
  for (uint32_t k = 1U; k < real_spectrum_size; k++) {
    max = fmaxf(spectrum[k], max);
  }
  return max;
}

float min_spectral_value(const float *spectrum, const uint32_t real_spectrum_size) {
  if (!spectrum || real_spectrum_size <= 0U) {
    return 0.F;
  }

  float min = spectrum[0];
  for (uint32_t k = 1U; k < real_spectrum_size; k++) {
    min = fminf(spectrum[k], min);
  }
  return min;
}

bool direct_matrix_to_vector_spectral_convolution(const float *matrix_spectum,
                                                  const float *spectrum,
                                                  float *out_spectrum,
                                                  uint32_t spectrum_size) {
  if (!matrix_spectum || !spectrum || !out_spectrum || spectrum_size <= 0) {
    return false;
  }

  for (uint32_t i = 0U; i < spectrum_size; i++) {
    out_spectrum[i] = 0.F;
    for (uint32_t j = 0U; j < spectrum_size; j++) {
      out_spectrum[i] += matrix_spectum[i * spectrum_size + j] * spectrum[j];
    }
  }

  return true;
}

inline float fft_bin_to_freq(const uint32_t bin_index,
                             const uint32_t sample_rate,
                             const uint32_t fft_size) {
  return (float)bin_index * ((float)sample_rate / (float)fft_size);
}

inline uint32_t freq_to_fft_bin(const float freq, const uint32_t sample_rate,
                                const uint32_t fft_size) {
  return (uint32_t)(freq / ((float)sample_rate / (float)fft_size / 2.F));
}

float spectral_flux(const float *spectrum, const float *previous_spectrum,
                    const uint32_t spectrum_size) {
  if (!spectrum || !previous_spectrum || spectrum_size <= 0U) {
    return 0.F;
  }

  float spectral_flux = 0.F;

  for (uint32_t i = 0U; i < spectrum_size; i++) {
    const float temp = sqrtf(spectrum[i]) - sqrtf(previous_spectrum[i]);
    spectral_flux += (temp + fabsf(temp)) / 2.F;
  }
  return spectral_flux;
}

bool get_rolling_mean_spectrum(float *averaged_spectrum,
                               const float *current_spectrum,
                               const uint32_t number_of_blocks,
                               const uint32_t spectrum_size) {
  if (!averaged_spectrum || !current_spectrum || (spectrum_size <= 0U)) {
    return false;
  }

  for (uint32_t k = 1U; k < spectrum_size; k++) {
    if (number_of_blocks <= 1U) {
      averaged_spectrum[k] = current_spectrum[k];
    } else {
      averaged_spectrum[k] += (averaged_spectrum[k] - current_spectrum[k]) /
                              (float)number_of_blocks;
    }
  }

  return true;
}