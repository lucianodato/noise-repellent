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

#ifndef NOISE_ESTIMATOR_H
#define NOISE_ESTIMATOR_H

#include <stdbool.h>

typedef struct NoiseEstimator NoiseEstimator;

bool is_noise_estimation_available(NoiseEstimator *self);
float *noise_estimation_run(NoiseEstimator *self, float *spectrum);
void noise_estimation_reset(NoiseEstimator *self);
NoiseEstimator *noise_estimation_initialize(int fft_size);
void noise_estimation_free(NoiseEstimator *self);

#endif