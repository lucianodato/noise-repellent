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

#ifndef SPL_SPECTRUM_CONVERTER_H
#define SPL_SPECTRUM_CONVERTER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct SplSpectrumConverter SplSpectrumConverter;

SplSpectrumConverter *reference_spectrum_initialize(uint32_t fft_size,
                                                    uint32_t sample_rate);
void reference_spectrum_free(SplSpectrumConverter *self);
bool convert_spectrum_to_dbspl(SplSpectrumConverter *self, float *spectrum);

#endif