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

#ifndef SPECTRAL_FEATURES_H
#define SPECTRAL_FEATURES_H

#include <stdbool.h>
#include <stdint.h>

typedef struct SpectralFeatures SpectralFeatures;

struct SpectralFeatures {
  float *power_spectrum;
  float *phase_spectrum;
  float *magnitude_spectrum;
};

bool get_fft_power_spectrum(const float *fft_spectrum,
                            uint32_t fft_spectrum_size,
                            float *fft_power_spectrum,
                            uint32_t fft_power_spectrum_size);
bool get_fft_magnitude_spectrum(const float *fft_spectrum,
                                uint32_t fft_spectrum_size,
                                float *fft_magnitude_spectrum,
                                uint32_t fft_magnitude_spectrum_size);
bool get_fft_phase_spectrum(const float *fft_spectrum,
                            uint32_t fft_spectrum_size,
                            float *fft_phase_spectrum,
                            uint32_t fft_phase_spectrum_size);

#endif