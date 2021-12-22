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

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Stft configurations

// Fft configurations
#define FFT_SIZE 2048

// OverlapAdd configurations
#define OVERLAP_FACTOR 2

// Windows
#define INPUT_WINDOW_TYPE 3
#define OUTPUT_WINDOW_TYPE 3

// ------------------------

// Denoiser configurations

// Transient protection
#define UPPER_LIMIT 5.f

// dB to dBSPL converter
#define REFERENCE_SINE_WAVE_FREQ 1000.f
#define REFERENCE_LEVEL 90.f
#define SINE_AMPLITUDE 1.f

// Spectral Whitening
#define WHITENING_DECAY_RATE 1000.f
#define WHITENING_FLOOR 0.02f

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

// ------------------------

#endif
