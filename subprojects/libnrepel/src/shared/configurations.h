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

#ifndef MODULES_CONFIGURATIONS_H
#define MODULES_CONFIGURATIONS_H

#include "critical_bands.h"
#include "fft_transform.h"
#include "spectral_features.h"
#include "spectral_utils.h"

#ifndef M_PI
#define M_PI 3.1415926535F
#endif

/* --------------------------------------------------------------------- */
/* ------------------- Shared Modules configurations ------------------- */
/* --------------------------------------------------------------------- */

// FFT configurations
#define PADDING_CONFIGURATION_GENERAL NO_PADDING
#define PADDING_CONFIGURATION_SPEECH FIXED_AMOUNT
#define ZEROPADDING_AMOUNT 25

// Absolute hearing thresholds
#define REFERENCE_SINE_WAVE_FREQ 1000.F
#define REFERENCE_LEVEL 90.F
#define SINE_AMPLITUDE 1.F

// Spectral Whitening
#define WHITENING_DECAY_RATE 1000.F
#define WHITENING_FLOOR 0.02F

// Masking Thresholds
#define BIAS 0
#define HIGH_FREQ_BIAS 20.F
#if BIAS
#define relative_thresholds                                                    \
  [N_BARK_BANDS] = {-16.F, -17.F, -18.F, -19.F, -20.F, -21.F, -22.F,           \
                    -23.F, -24.F, -25.F, -25.F, -25.F, -25.F, -25.F,           \
                    -25.F, -24.F, -23.F, -22.F, -19.F, -18.F, -18.F,           \
                    -18.F, -18.F, -18.F, -18.F}
#endif

// Gain Estimators
#define GAMMA1 2.F
#define GAMMA2 0.5F

#define ALPHA_MAX 6.F
#define ALPHA_MIN 1.F
#define BETA_MAX 0.02F
#define BETA_MIN 0.F

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

/* ----------------------------------------------------------- */
/* ------------------- Stft configurations ------------------- */
/* ----------------------------------------------------------- */

// STFT configurations - Frame size in milliseconds
#define FRAME_SIZE_GENERAL 46
#define FRAME_SIZE_SPEECH 20

// OverlapAdd configurations
#define OVERLAP_FACTOR_GENERAL 4
#define OVERLAP_FACTOR_SPEECH 2

// Windows
#define INPUT_WINDOW_TYPE_GENERAL VORBIS_WINDOW
#define OUTPUT_WINDOW_TYPE_GENERAL VORBIS_WINDOW

#define INPUT_WINDOW_TYPE_SPEECH VORBIS_WINDOW
#define OUTPUT_WINDOW_TYPE_SPEECH VORBIS_WINDOW

/* --------------------------------------------------------------- */
/* ------------------- Denoiser configurations ------------------- */
/* --------------------------------------------------------------- */

// Spectral Type
#define SPECTRAL_TYPE_GENERAL POWER_SPECTRUM

// Transient protection
#define UPPER_LIMIT 5.F

// Masking
#define N_CRITICAL_BANDS 25
#define CRITICAL_BANDS_TYPE BARK_SCALE

/* ------------------------------------------------------------------------ */
/* ------------------- Adaptive Denoiser configurations ------------------- */
/* ------------------------------------------------------------------------ */

// Spectral Type
#define SPECTRAL_TYPE_SPEECH POWER_SPECTRUM

// Masking
#define N_CRITICAL_BANDS_SPEECH 25
#define CRITICAL_BANDS_TYPE_SPEECH BARK_SCALE
#define DEFAULT_MASKING_CEILING 2.F
#define DEFAULT_MASKING_FLOOR 0.01F

#endif // ifndef
