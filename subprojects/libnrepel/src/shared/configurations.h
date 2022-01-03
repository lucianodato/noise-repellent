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

#include "spectral_features.h"
#include "spectral_utils.h"

#ifndef M_PI
#define M_PI 3.1415926535F
#endif

/* --------------------------------------------------------------------- */
/* ------------------- Shared Modules configurations ------------------- */
/* --------------------------------------------------------------------- */

// Fft configurations
#define FFT_SIZE 2048

// dB to dBSPL converter
#define REFERENCE_SINE_WAVE_FREQ 1000.F
#define REFERENCE_LEVEL 90.F
#define SINE_AMPLITUDE 1.F

// Spectral Whitening
#define WHITENING_DECAY_RATE 1000.F
#define WHITENING_FLOOR 0.02F

/* ----------------------------------------------------------- */
/* ------------------- Stft configurations ------------------- */
/* ----------------------------------------------------------- */

// OverlapAdd configurations
#define OVERLAP_FACTOR 2U

// Windows
#define INPUT_WINDOW_TYPE VORBIS_WINDOW
#define OUTPUT_WINDOW_TYPE VORBIS_WINDOW

/* --------------------------------------------------------------- */
/* ------------------- Denoiser configurations ------------------- */
/* --------------------------------------------------------------- */

// Spectral Denoiser
#define SPECTRAL_TYPE POWER_SPECTRUM

// Transient protection
#define UPPER_LIMIT 5.F

// Masking
#define N_BARK_BANDS 25U
#define BIAS 0U
#define HIGH_FREQ_BIAS 20.F

#if BIAS
#define relative_thresholds                                                    \
  [N_BARK_BANDS] = {-16.F, -17.F, -18.F, -19.F, -20.F, -21.F, -22.F,           \
                    -23.F, -24.F, -25.F, -25.F, -25.F, -25.F, -25.F,           \
                    -25.F, -24.F, -23.F, -22.F, -19.F, -18.F, -18.F,           \
                    -18.F, -18.F, -18.F, -18.F}
#endif

// Gain Estimator
#define GAMMA1 2.F
#define GAMMA2 0.5F

#define ALPHA_MAX 6.F
#define ALPHA_MIN 1.F
#define BETA_MAX 0.02F
#define BETA_MIN 0.F

// Noise Estimator
#define MIN_NUMBER_OF_WINDOWS_NOISE_AVERAGED 5U

// Louizou Estimator
#define N_SMOOTH 0.7F
#define BETA_AT 0.8F
#define GAMMA 0.998F
#define ALPHA_P 0.2F
#define ALPHA_D 0.95F

#define CROSSOVER_POINT1 1000.F
#define CROSSOVER_POINT2 3000.F
#define BAND_1_GAIN 2.F
#define BAND_2_GAIN 2.F
#define BAND_3_GAIN 7.F

#endif // ifndef
