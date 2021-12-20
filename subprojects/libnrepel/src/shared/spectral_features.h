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

SpectralFeatures *spectral_features_initialize(uint32_t spectral_size);
void spectral_features_free(SpectralFeatures *self);

bool compute_power_spectrum(SpectralFeatures *self, const float *fft_spectrum,
                            uint32_t fft_spectrum_size);
bool compute_magnitude_spectrum(SpectralFeatures *self,
                                const float *fft_spectrum,
                                uint32_t fft_spectrum_size);
bool compute_phase_spectrum(SpectralFeatures *self, const float *fft_spectrum,
                            uint32_t fft_spectrum_size);
float *get_power_spectrum(SpectralFeatures *self);
float *get_magnitude_spectrum(SpectralFeatures *self);
float *get_phase_spectrum(SpectralFeatures *self);

#endif