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

#ifndef MODULES_CONFIGURATIONS_H
#define MODULES_CONFIGURATIONS_H

#include "critical_bands.h"
#include "fft_transform.h"
#include "gain_estimators.h"
#include "noise_scaling_criterias.h"
#include "spectral_features.h"
#include "spectral_utils.h"
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.1415926535F
#endif

/* --------------------------------------------------------------------- */
/* ------------------- Shared Modules configurations ------------------- */
/* --------------------------------------------------------------------- */

// FFT configurations
#define ZEROPADDING_AMOUNT 50 // Even Number

// Absolute hearing thresholds
#define REFERENCE_SINE_WAVE_FREQ 1000.F
#define REFERENCE_LEVEL 90.F
#define SINE_AMPLITUDE 1.F

// Spectral Whitening
#define WHITENING_DECAY_RATE 1000.F
#define WHITENING_FLOOR 0.02F

// Masking Thresholds
#define BIAS false
#define HIGH_FREQ_BIAS 20.F
#if BIAS
#define relative_thresholds                                                    \
  [N_BARK_BANDS] = {-16.F, -17.F, -18.F, -19.F, -20.F, -21.F, -22.F,           \
                    -23.F, -24.F, -25.F, -25.F, -25.F, -25.F, -25.F,           \
                    -25.F, -24.F, -23.F, -22.F, -19.F, -18.F, -18.F,           \
                    -18.F, -18.F, -18.F, -18.F}
#endif

// Postfilter SNR Threshold
#define POSTFILTER_THRESHOLD 0.4F
#define POSTFILTER_SCALE 10.0F
#define PRESERVE_MINIMUN_GAIN true

// Gain Estimators
#define GSS_EXPONENT                                                           \
  2.0F // 2 Power Subtraction / 1 Magnitude Subtraxtion / 0.5 Spectral
       // Subtraction

// Oversubtraction criterias
#define ALPHA_MAX 6.F
#define ALPHA_MIN 1.F
#define BETA_MAX 0.02F
#define BETA_MIN 0.F
#define DEFAULT_OVERSUBTRACTION ALPHA_MIN
#define DEFAULT_UNDERSUBTRACTION BETA_MAX

// Noise Estimator
#define MIN_NUMBER_OF_WINDOWS_NOISE_AVERAGED 5

// Adaptive Estimator
#define N_SMOOTH 0.7F
#define BETA_AT 0.8F
#define GAMMA 0.998F
#define ALPHA_P 0.2F
#define ALPHA_D 0.85F

#define CROSSOVER_POINT1 1000.F
#define CROSSOVER_POINT2 3000.F
#define BAND_1_LEVEL 2.F
#define BAND_2_LEVEL 2.F
#define BAND_3_LEVEL 5.F

/* --------------------------------------------------------------- */
/* ------------------- Denoiser configurations ------------------- */
/* --------------------------------------------------------------- */

// STFT configurations - Frame size in milliseconds
#define FRAME_SIZE_GENERAL 46
#define OVERLAP_FACTOR_GENERAL 4
#define INPUT_WINDOW_TYPE_GENERAL VORBIS_WINDOW
#define OUTPUT_WINDOW_TYPE_GENERAL VORBIS_WINDOW

// Fft configuration
#define PADDING_CONFIGURATION_GENERAL NO_PADDING

// Spectral Type
#define SPECTRAL_TYPE_GENERAL POWER_SPECTRUM

// Transient protection
#define UPPER_LIMIT 5.F

// Masking
#define CRITICAL_BANDS_TYPE BARK_SCALE

// Noise Scaling strategy
#define OVERSUBTRACTION_TYPE MASKING_THRESHOLDS
#define GAIN_ESTIMATION_TYPE WIENER

/* ------------------------------------------------------------------------ */
/* ------------------- Adaptive Denoiser configurations ------------------- */
/* ------------------------------------------------------------------------ */

// STFT configurations - Frame size in milliseconds
#define FRAME_SIZE_SPEECH 26
#define OVERLAP_FACTOR_SPEECH 2
#define INPUT_WINDOW_TYPE_SPEECH VORBIS_WINDOW
#define OUTPUT_WINDOW_TYPE_SPEECH VORBIS_WINDOW

// Spectral Type
#define SPECTRAL_TYPE_SPEECH POWER_SPECTRUM

// Fft configurations
#define PADDING_CONFIGURATION_SPEECH NO_PADDING

// Masking
#define CRITICAL_BANDS_TYPE_SPEECH OPUS_SCALE

// Noise Scaling strategy
#define OVERSUBTRACTION_TYPE_SPEECH MASKING_THRESHOLDS
#define GAIN_ESTIMATION_TYPE_SPEECH WIENER

#endif // ifndef
