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

float bin_to_freq(int i, float samp_rate, int N);
void compute_bark_mapping(MaskingEstimator *self);
void compute_absolute_thresholds(MaskingEstimator *self);
void hanning_window(float *window, int N);
void get_power_spectrum(MaskingEstimator *self, float *window, float *signal, float *power_spectrum);
void spl_reference(MaskingEstimator *self);
void compute_spectral_spreading_function(MaskingEstimator *self);
void convolve_with_spectral_spreading_function(MaskingEstimator *self, float *bark_spectrum, float *spreaded_spectrum);
void compute_bark_spectrum(MaskingEstimator *self, float *bark_spectrum, float *spectrum,
						   float *intermediate_band_bins, float *n_bins_per_band);
void convert_to_dbspl(MaskingEstimator *self, float *masking_thresholds);
float compute_tonality_factor(float *spectrum, float *intermediate_band_bins,
							  float *n_bins_per_band, int band);
void compute_masking_thresholds(MaskingEstimator *self, float *spectrum, float *masking_thresholds);
void masking_estimation_reset(MaskingEstimator *self);
MaskingEstimator *masking_estimation_initialize(int fft_size, int samp_rate);
void masking_estimation_free(MaskingEstimator *self);

#endif