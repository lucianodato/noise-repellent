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

#ifndef SPECTRAL_PROCESSOR_H
#define SPECTRAL_PROCESSOR_H

#include <stdbool.h>

// Generic Spectral Processing function over an FFT spectrum. Receives any
// spectral processing module handle (void *) and the FFT of a audio block.
// This is to inject any spectral processor and processing function into the
// STFT transform at runtime
typedef void *SpectralProcessorHandle;

// Processing function which deals with the fft spectrum by mutating the array
// with any DSP that operates with the FFT spectrum (1d FFTW spectrum)
typedef bool (*spectral_processing)(SpectralProcessorHandle spectral_processor,
                                    float *fft_spectrum);
#endif