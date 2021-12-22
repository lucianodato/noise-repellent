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

#ifndef STFT_PROCESSOR_H
#define STFT_PROCESSOR_H

#include <stdbool.h>
#include <stdint.h>

typedef struct StftProcessor StftProcessor;

// Generic Spectral Processing function over an FFT spectrum. Receives any
// spectral processing module handle (void *) and the FFT of a audio block.
typedef bool spectral_processing(void *spectral_processor, float *fft_spectrum);

StftProcessor *stft_processor_initialize();
void stft_processor_free(StftProcessor *self);
uint32_t get_stft_latency(StftProcessor *self);
bool stft_processor_run(StftProcessor *self,
                        spectral_processing *spectral_processing,
                        void *spectral_processor, uint32_t number_of_samples,
                        const float *input, float *output);
uint32_t get_fft_size(StftProcessor *self);
uint32_t get_overlap_factor(StftProcessor *self);
uint32_t get_spectral_processing_size(StftProcessor *self);

#endif