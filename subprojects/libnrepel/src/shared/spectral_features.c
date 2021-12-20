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

#include "spectral_features.h"
#include <math.h>

bool get_fft_power_spectrum(const float *fft_spectrum,
                            const uint32_t fft_spectrum_size,
                            float *fft_power_spectrum,
                            const uint32_t fft_power_spectrum_size) {
  if (!fft_spectrum || !fft_power_spectrum || !fft_spectrum_size ||
      !fft_power_spectrum_size) {
    return false;
  }

  float real_bin = fft_spectrum[0];

  fft_power_spectrum[0] = real_bin * real_bin;

  for (uint32_t k = 1; k <= fft_power_spectrum_size; k++) {
    float power = 0.f;

    real_bin = fft_spectrum[k];
    float imag_bin = fft_spectrum[fft_spectrum_size - k];

    if (k < fft_power_spectrum_size) {
      power = (real_bin * real_bin + imag_bin * imag_bin);
    } else {
      power = real_bin * real_bin;
    }

    fft_power_spectrum[k] = power;
  }

  return true;
}

bool get_fft_magnitude_spectrum(const float *fft_spectrum,
                                const uint32_t fft_spectrum_size,
                                float *fft_magnitude_spectrum,
                                const uint32_t fft_magnitude_spectrum_size) {

  if (!fft_spectrum || !fft_magnitude_spectrum || !fft_spectrum_size ||
      !fft_magnitude_spectrum_size) {
    return false;
  }

  float real_bin = fft_spectrum[0];

  fft_magnitude_spectrum[0] = real_bin;

  for (uint32_t k = 1; k <= fft_magnitude_spectrum_size; k++) {
    float magnitude = 0.f;

    real_bin = fft_spectrum[k];
    float imag_bin = fft_spectrum[fft_spectrum_size - k];

    if (k < fft_magnitude_spectrum_size) {
      magnitude = sqrtf(real_bin * real_bin + imag_bin * imag_bin);

    } else {
      magnitude = real_bin;
    }

    fft_magnitude_spectrum[k] = magnitude;
  }

  return true;
}

bool get_fft_phase_spectrum(const float *fft_spectrum,
                            uint32_t fft_spectrum_size,
                            float *fft_phase_spectrum,
                            uint32_t fft_phase_spectrum_size) {
  if (!fft_spectrum || !fft_phase_spectrum || !fft_spectrum_size ||
      !fft_phase_spectrum_size) {
    return false;
  }

  float real_bin = fft_spectrum[0];
  fft_phase_spectrum[0] = atan2f(real_bin, 0.f);

  for (uint32_t k = 1; k <= fft_phase_spectrum_size; k++) {
    float phase = 0.f;

    real_bin = fft_spectrum[k];
    float imag_bin = fft_spectrum[fft_spectrum_size - k];

    if (k < fft_phase_spectrum_size) {
      phase = atan2f(real_bin, imag_bin);
    } else {
      phase = atan2f(real_bin, 0.f);
    }

    fft_phase_spectrum[k] = phase;
  }

  return true;
}