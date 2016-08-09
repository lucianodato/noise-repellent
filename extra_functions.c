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

inline float Index2Freq(int i, float samp_rate, int N) {
  return (float) i * (samp_rate / N / 2.f);
}

inline int Freq2Index(float freq, float samp_rate, int N) {
  return (int) (freq / (samp_rate / N / 2.f));
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

inline float blackman(int k, int N) {
  float p = ((float)(k))/((float)(N-1));
  return 0.42-0.5*cosf(2.f*M_PI*p) + 0.08*cosf(4.f*M_PI*p);
}

inline float hanning(int k, int N) {
  float p = ((float)(k))/((float)(N-1));
  return 0.5 - 0.5 * cosf(2.f*M_PI*p);
}

inline float hamming(int k, int N) {
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

static float max_spectral_value(float* noise_print, int N){
  int k;
  float max = 0.f;
  for(k = 0; k <= N; k++){
    max = MAX(noise_print[k],max);
  }
  return max;
}

static float min_spectral_value(float* noise_print, int N){
  int k;
  float min = FLT_MAX;
  for(k = 0; k <= N; k++){
    min = MIN(noise_print[k],min);
  }
  return min;
}

//unnormalized Hann windows for whitening tappering
void tappering_filter_calc(float* filter, int N,float WA) {
  int k;
  for (k = 0; k < N; k++){
    filter[k] = powf(hanning(k, N),WA);//Half hann window tappering in favor of high frequencies
  }
}

void whitening_of_spectrum(float* spectrum,float wa,int N){
  float whitened_spectrum[N];
  float tmp_min = min_spectral_value(spectrum,N);
  float tmp_max = max_spectral_value(spectrum,N);
  for (int k = 0; k <= N; k++) {
    if(spectrum[k] > FLT_MIN){ //Protects against division by 0
      whitened_spectrum[k] = powf((spectrum[k]/(tmp_max-tmp_min)),wa);
      spectrum[k] /= whitened_spectrum[k];
      if(k < N){
        spectrum[N-k] /= whitened_spectrum[k];
      }
    }
  }
}

void apply_tappering_filter(float* spectrum,float* filter,int N) {
  for (int k = 0; k <= N; k++) {
    if(spectrum[k] > FLT_MIN) {
      spectrum[k] *= filter[k];//Half hann window tappering in favor of high frequencies
      if(k < N) {
        spectrum[N-k] *= filter[k];//Half hann window tappering in favor of high frequencies
      }
    }
  }
}

// //Spectral smoothing (Based on Audacity code)
void spectral_smoothing(float* spectrum, int smoothing_bins,int N, int middle_bin){
  int k;
  float smoothing_tmp[N+1];
  float log_spectrum[N+1];

  //Initialize smothingbins_tmp
  for (k = 0; k <= N; ++k) {
    log_spectrum[k] = logf(spectrum[k]);
    smoothing_tmp[k] = 0.f;
  }

  //do not smooth up to the middle bin not even DC
  for (k = 0; k < middle_bin; ++k) {
    smoothing_tmp[k] = log_spectrum[k];
  }

  for (k = middle_bin; k < N; ++k) {
    const int j0 = MAX(middle_bin, k - smoothing_bins);
    const int j1 = MIN(N, k + smoothing_bins);
    for(int l = j0; l <= j1; ++l) {
       smoothing_tmp[k] += log_spectrum[l];
    }
    log_spectrum[k] /= (j1 - j0 + 1);
  }

  for (k = 0; k <= N; ++k){
      spectrum[k] = expf(smoothing_tmp[k]);
  }
}
