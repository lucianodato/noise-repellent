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

#ifndef DENOISE_MIXER_H
#define DENOISE_MIXER_H

#include <stdbool.h>
#include <stdint.h>

typedef struct DenoiseMixerParameters {
  float noise_level;
  bool residual_listen;
  float whitening_amount;
} DenoiseMixerParameters;

typedef struct DenoiseMixer DenoiseMixer;

DenoiseMixer *denoise_mixer_initialize(uint32_t fft_size, uint32_t sample_rate,
                                       uint32_t hop);
void denoise_mixer_free(DenoiseMixer *self);
bool denoise_mixer_run(DenoiseMixer *self, float *fft_spectrum,
                       const float *gain_spectrum,
                       DenoiseMixerParameters parameters);

#endif