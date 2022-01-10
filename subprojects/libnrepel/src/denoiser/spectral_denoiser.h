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

#ifndef SPECTRAL_DENOISER_H
#define SPECTRAL_DENOISER_H

#include "../../include/nrepel.h"
#include "../shared/noise_profile.h"
#include "../shared/spectral_processor.h"
#include <stdbool.h>
#include <stdint.h>

SpectralProcessorHandle spectral_denoiser_initialize(
    uint32_t sample_rate, uint32_t fft_size, uint32_t overlap_factor,
    NoiseProfile *noise_profile, NrepelDenoiseParameters *parameters);
void spectral_denoiser_free(SpectralProcessorHandle instance);
bool spectral_denoiser_run(SpectralProcessorHandle instance,
                           float *fft_spectrum);
bool spectral_adaptive_denoiser_run(SpectralProcessorHandle instance,
                                    float *fft_spectrum);

#endif