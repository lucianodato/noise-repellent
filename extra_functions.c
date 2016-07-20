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


#include <math.h>
#include <cmath>
#include <float.h>


//Window types
#define HANNING_WINDOW 0
#define HAMMING_WINDOW 1

//AUXILIARY Functions

// Force already-denormal float value to zero
static inline float sanitize_denormal(float value) {
    if (isnan(value)) {
      return FLT_MIN;
      //return 0.f;
    } else {
      return value;
    }

}

static inline int sign(float x) {
        return (x >= 0.f ? 1 : -1);
}

static inline float from_dB(float gdb) {
        return (exp(gdb/20.f*log(10.f)));
}

static inline float to_dB(float g) {
        return (20.f*log10(g));
}

static float hanning(int k, int N) {
  float p = ((float)(k))/((float)(N));
  return 0.5 - 0.5 * cos(2.f*M_PI*p);
}

static float hamming(int k, int N) {
  float p = ((float)(k))/((float)(N));
  return 0.54 - 0.46 * cos(2.f*M_PI*p);
}

static void fft_window(float* window, int N, int window_type) {
  float value = 0.f;
  float sum_values = 0.f;
  int k;
  for (k = 0; k < N; k++){
    switch (window_type){
      case HANNING_WINDOW:
        value = hanning(k, N);
      break;
      case HAMMING_WINDOW:
        value = hamming(k, N);
      break;
    }
    window[k] = value;
    sum_values += value;
  }

  for (k = 0; k < N; k++){
    window[k] /= sum_values; //Normalized Window
  }
}
