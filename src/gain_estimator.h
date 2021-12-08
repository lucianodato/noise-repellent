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

#ifndef GAIN_ESTIMATOR_H
#define GAIN_ESTIMATOR_H

typedef struct GainEstimator GainEstimator;

float max_spectral_value(float *spectrum, int N);
float min_spectral_value(float *spectrum, int N);
void wiener_subtraction(GainEstimator *self);
void spectral_gating(GainEstimator *self);
void denoise_gain_generalized_spectral_substraction(GainEstimator *self);
void compute_alpha_and_beta(GainEstimator *self, float masking_ceiling_limit, float masking_floor_limit);
void gain_estimation_run(GainEstimator *self, float *signal_spectrum, float *noise_profile, float *gain_spectrum, float transient_threshold,
						 float masking_ceiling_limit, float release, float noise_rescale);
void gain_estimation_free(GainEstimator *self);
GainEstimator *gain_estimation_initialize(int fft_size, int samp_rate, int hop);

#endif