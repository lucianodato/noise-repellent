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

/**
* \file spectrum_smoother.c
* \author Luciano Dato
* \brief Contains a spectrum smoother abstraction
*/

#ifndef SPECTRUM_SMOOTHER_H
#define SPECTRUM_SMOOTHER_H

typedef struct SpectralSmoother SpectralSmoother;

void get_release_coefficient(SpectralSmoother *self, float release);
void apply_time_envelope(SpectralSmoother *self);
void spectral_smoothing_run(SpectralSmoother *self, float release);
void spectral_smoothing_reset(SpectralSmoother *self);
SpectralSmoother *spectral_smoothing_initialize(int fft_size, int samp_rate, int hop);
void spectral_smoothing_free(SpectralSmoother *self);

#endif