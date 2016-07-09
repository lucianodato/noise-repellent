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
#include <assert.h>

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

//This function was done by Ethan Brodsky
void unwrap(float* p, int N)
 // ported from matlab (Dec 2002)
  {
    float dp[MAX_LENGTH];
    float dps[MAX_LENGTH];
    float dp_corr[MAX_LENGTH];
    float cumsum[MAX_LENGTH];
    float cutoff = M_PI;               /* default value in matlab */
    int j;

    assert(N <= MAX_LENGTH);

    //Fill dp_corr with zeros to avoid warnings from the compiler
    for (j = 0; j<MAX_LENGTH;j++){
      dp_corr[j]=0;
    }

   // incremental phase variation
   // MATLAB: dp = diff(p, 1, 1);
    for (j = 0; j < N-1; j++)
      dp[j] = p[j+1] - p[j];

   // equivalent phase variation in [-pi, pi]
   // MATLAB: dps = mod(dp+dp,2*pi) - pi;
    for (j = 0; j < N-1; j++)
      dps[j] = (dp[j]+M_PI) - floor((dp[j]+M_PI) / (2*M_PI))*(2*M_PI) - M_PI;

   // preserve variation sign for +pi vs. -pi
   // MATLAB: dps(dps==pi & dp>0,:) = pi;
    for (j = 0; j < N-1; j++)
      if ((dps[j] == -M_PI) && (dp[j] > 0))
        dps[j] = M_PI;

   // incremental phase correction
   // MATLAB: dp_corr = dps - dp;
    for (j = 0; j < N-1; j++)
      dp_corr[j] = dps[j] - dp[j];

   // Ignore correction when incremental variation is smaller than cutoff
   // MATLAB: dp_corr(abs(dp)<cutoff,:) = 0;
    for (j = 0; j < N-1; j++)
      if (fabs(dp[j]) < cutoff)
        dp_corr[j] = 0;

   // Find cumulative sum of deltas
   // MATLAB: cumsum = cumsum(dp_corr, 1);
    cumsum[0] = dp_corr[0];
    for (j = 1; j < N-1; j++)
      cumsum[j] = cumsum[j-1] + dp_corr[j];

   // Integrate corrections and add to P to produce smoothed phase values
   // MATLAB: p(2:m,:) = p(2:m,:) + cumsum(dp_corr,1);
    for (j = 1; j < N; j++)
      p[j] += cumsum[j-1];
  }
