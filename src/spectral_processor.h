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

#ifndef SPECTRAL_PROCESSOR_H
#define SPECTRAL_PROCESSOR_H

#include "internal_data_types.h"
#include <float.h>
#include <stdbool.h>
#include <stdint.h>

typedef struct SpectralProcessor SpectralProcessor;

SpectralProcessor *spectral_processor_initialize(uint32_t sample_rate,
                                                 uint32_t fft_size,
                                                 uint32_t overlap_factor);
void spectral_processor_free(SpectralProcessor *self);
void load_processor_parameters(SpectralProcessor *self,
                               ProcessorParameters *new_parameters);
void spectral_processor_run(SpectralProcessor *self, float *fft_spectrum);

// Polymorphysm with this one
void load_noise_profile(SpectralProcessor *self, NoiseProfile *noise_profile);

#endif