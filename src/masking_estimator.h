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

#ifndef MASKING_ESTIMATOR_H
#define MASKING_ESTIMATOR_H

typedef struct MaskingEstimator MaskingEstimator;

MaskingEstimator *masking_estimation_initialize(int fft_size, int samp_rate);
void masking_estimation_free(MaskingEstimator *self);
void compute_masking_thresholds(MaskingEstimator *self, float *spectrum, float *masking_thresholds);
void compute_bark_mapping(MaskingEstimator *self);
void compute_absolute_thresholds(MaskingEstimator *self);
void spl_reference(MaskingEstimator *self);
void compute_spectral_spreading_function(MaskingEstimator *self);
void convolve_with_spectral_spreading_function(MaskingEstimator *self, float *bark_spectrum, float *spreaded_spectrum);

#endif