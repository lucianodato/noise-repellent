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

#ifndef GAIN_ESTIMATORS_H
#define GAIN_ESTIMATORS_H

#include <stdbool.h>
#include <stdint.h>

void denoise_mixer(uint32_t fft_size, float *fft_spectrum,
                   const float *gain_spectrum, float *denoised_spectrum,
                   float *residual_spectrum, bool residual_listen,
                   float reduction_amount);
void estimate_gains(uint32_t real_spectrum_size, uint32_t fft_size,
                    const float *spectrum, const float *noise_spectrum,
                    float *gain_spectrum, const float *alpha,
                    const float *beta);
void spectral_gating(uint32_t real_spectrum_size, uint32_t fft_size,
                     const float *spectrum, float *noise_spectrum,
                     float *gain_spectrum, const float *alpha);
void wiener_subtraction(uint32_t real_spectrum_size, uint32_t fft_size,
                        const float *spectrum, float *noise_spectrum,
                        float *gain_spectrum, const float *alpha);
void generalized_spectral_subtraction(uint32_t real_spectrum_size,
                                      uint32_t fft_size, const float *spectrum,
                                      const float *noise_spectrum,
                                      float *gain_spectrum, const float *alpha,
                                      const float *beta);

#endif