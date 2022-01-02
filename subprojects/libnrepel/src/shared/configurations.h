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
#define M_PI (3.14159265358979323846F)
#endif

/* --------------------------------------------------------------------- */
/* ------------------- Shared Modules configurations ------------------- */
/* --------------------------------------------------------------------- */

// Fft configurations
#define FFT_SIZE 2048

// dB to dBSPL converter
#define REFERENCE_SINE_WAVE_FREQ 1000.f
#define REFERENCE_LEVEL 90.f
#define SINE_AMPLITUDE 1.f

// Spectral Whitening
#define WHITENING_DECAY_RATE 1000.f
#define WHITENING_FLOOR 0.02f

/* ----------------------------------------------------------- */
/* ------------------- Stft configurations ------------------- */
/* ----------------------------------------------------------- */

// OverlapAdd configurations
#define OVERLAP_FACTOR 2

// Windows
#define INPUT_WINDOW_TYPE VORBIS_WINDOW
#define OUTPUT_WINDOW_TYPE VORBIS_WINDOW

/* --------------------------------------------------------------- */
/* ------------------- Denoiser configurations ------------------- */
/* --------------------------------------------------------------- */

// Spectral Denoiser
#define SPECTRAL_TYPE POWER_SPECTRUM

// Transient protection
#define UPPER_LIMIT 5.f

// Masking
#define N_BARK_BANDS 25
#define BIAS 0
#define HIGH_FREQ_BIAS 20.f

#if BIAS
#define relative_thresholds                                                    \
  [N_BARK_BANDS] = {-16.f, -17.f, -18.f, -19.f, -20.f, -21.f, -22.f,           \
                    -23.f, -24.f, -25.f, -25.f, -25.f, -25.f, -25.f,           \
                    -25.f, -24.f, -23.f, -22.f, -19.f, -18.f, -18.f,           \
                    -18.f, -18.f, -18.f, -18.f}
#endif

// Gain Estimator
#define GAMMA1 2.f
#define GAMMA2 0.5f

#define ALPHA_MAX 6.f
#define ALPHA_MIN 1.f
#define BETA_MAX 0.02f
#define BETA_MIN 0.0f

// Noise Estimator
#define MIN_NUMBER_OF_WINDOWS_NOISE_AVERAGED 5

// Louizou Estimator
#define N_SMOOTH 0.7F
#define BETA_AT 0.8F
#define GAMMA 0.998F
#define ALPHA_P 0.2F
#define ALPHA_D 0.95F

#define CROSSOVER_POINT1 1000.F
#define CROSSOVER_POINT2 3000.F
#define BAND_1_GAIN 2.0F
#define BAND_2_GAIN 2.0F
#define BAND_3_GAIN 7.0F

#endif // ifndef
