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

#ifndef GAIN_ESTIMATOR_H
#define GAIN_ESTIMATOR_H

#include "../../../include/nrepel.h"
#include "../../shared/noise_profile.h"
#include <stdbool.h>
#include <stdint.h>

typedef struct GainEstimator GainEstimator;

GainEstimator *gain_estimation_initialize(uint32_t fft_size,
                                          uint32_t sample_rate, uint32_t hop,
                                          NrepelDenoiseParameters *parameters,
                                          NoiseProfile *noise_profile);
void gain_estimation_free(GainEstimator *self);
bool gain_estimation_run(GainEstimator *self, const float *signal_spectrum,
                         float *gain_spectrum);
bool gain_estimation_run_adaptive(GainEstimator *self,
                                  const float *signal_spectrum,
                                  float *gain_spectrum);

#endif