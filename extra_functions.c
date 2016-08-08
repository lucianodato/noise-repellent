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
#include <float.h>


//Window types
#define HANNING_WINDOW 0
#define HAMMING_WINDOW 1
#define BLACKMAN_WINDOW 2

//AUXILIARY Functions

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

// Force already-denormal float value to zero
inline float sanitize_denormal(float value) {
  if (isnan(value)) {
    return FLT_MIN; //to avoid log errors
    //return 0.f; //to avoid log errors
  } else {
    return value;
  }

}

inline int sign(float x) {
  return (x >= 0.f ? 1.f : -1.f);
}

inline float from_dB(float gdb) {
  return (expf(gdb/20.f*logf(10.f)));
}

inline float to_dB(float g) {
  return (20.f*log10f(g));
}

static float blackman(int k, int N) {
  float p = ((float)(k))/((float)(N-1));
  return 0.42-0.5*cosf(2.f*M_PI*p) + 0.08*cosf(4.f*M_PI*p);
}

static float hanning(int k, int N) {
  float p = ((float)(k))/((float)(N-1));
  return 0.5 - 0.5 * cosf(2.f*M_PI*p);
}

static float hamming(int k, int N) {
  float p = ((float)(k))/((float)(N-1));
  return 0.54 - 0.46 * cosf(2.f*M_PI*p);
}

void fft_window(float* window, int N, int window_type) {
  float sum_values = 0.f;
  int k;
  for (k = 0; k < N; k++){
    switch (window_type){
      case BLACKMAN_WINDOW:
      window[k] = blackman(k, N);
      break;
      case HANNING_WINDOW:
      window[k] = hanning(k, N);
      break;
      case HAMMING_WINDOW:
      window[k] = hamming(k, N);
      break;
    }
    sum_values += window[k];
  }

  for (k = 0; k < N; k++){
    window[k] /= sum_values; //Normalized Window
  }
}

//unnormalized Hann windows for whitening tappering
void tappering_filter_calc(float* filter, int N,float WA) {
  int k;
  for (k = 0; k < N; k++){
    filter[k] = powf(hanning(k, N),WA);//Half hann window tappering in favor of high frequencies
  }
}
