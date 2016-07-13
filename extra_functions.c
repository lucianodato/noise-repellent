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
#define BLACKMAN_WINDOW 2
#define MAX_LENGTH 10000 //Unwrap - depends on fft size

//AUXILIARY Functions

// Force already-denormal float value to zero
static inline void
sanitize_denormal(float& value) {
    if (isnan(value)) {
        value = 0.f;
    }
}

static inline int
sign(float x) {
        return (x >= 0.f ? 1 : -1);
}

static inline float
from_dB(float gdb) {
        return (exp(gdb/20.f*log(10.f)));
}

static inline float
to_dB(float g) {
        return (20.f*log10(g));
}

float blackman(int k, int N)
{
  float p = ((float)(k))/(float)(N-1);
  return 0.42-0.5*cos(2.0*M_PI*p) + 0.08*cos(4.0*M_PI*p);
}

float hanning(int k, int N)
{
  float p = ((float)(k))/(float)(N-1);
  return 0.5 - 0.5 * cos(2.0*M_PI*p);
}

float hamming(int k, int N) {
  float p = ((float)(k))/(float)(N-1);
  return 0.54 - (0.46 * cos(2.0*M_PI*p));
}

void fft_window(float* window,int N, int window_type)
{
  float value,sum;
  sum = 0;
  int k;
  for (k = 0; k < N; k++){
    switch (window_type){
      case BLACKMAN_WINDOW:
        value = blackman(k, N);
      break;
      case HANNING_WINDOW:
        value = hanning(k, N);
      break;
      case HAMMING_WINDOW:
        value = hamming(k, N);
      break;
      default:
        value = 0;
      }
      window[k]= value;
      sum += value;
  }
  //Normalize Window
  for (k = 0; k < N; k++){
    window[k]/=sum;
  }
}

float Atan2modif(float x, float y)
{
  float signx;
  if (x > 0.) signx = 1.;
  else signx = -1.;

  if (x == 0.) return 0.;
  if (y == 0.) return signx * M_PI / 2.;

  return atan2(x, y);
}
