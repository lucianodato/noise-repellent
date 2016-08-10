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

inline float mean(int m, float* a) {
    int sum=0, i;
    for(i=0; i<m; i++)
        sum+=a[i];
    return((float)sum/m);
}

inline float median(int n, float* x) {
    float temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.f);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

inline float moda(int n, float* x) {
  float temp[n];
  int i,j,pos_max;
  float max;

  for(i = 0;i<n; i++) {
      temp[i]=0.f;
  }

  for(i=0; i<n; i++) {
      for(j=i; j<n; j++) {
          if(x[j] == x[i]) temp[i]++;
      }
  }

  max=temp[0];
  pos_max = 0;
  for(i=0; i<n; i++) {
      if(temp[i] > max) {
          pos_max = i;
          max=temp[i];
      }
  }
  return x[pos_max];
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
    //   if(k < N){
    //     spectrum[N-k] /= whitened_spectrum[k];
    //   }
    }
  }
}

void apply_tappering_filter(float* spectrum,float* filter,int N) {
  for (int k = 0; k <= N; k++) {
    if(spectrum[k] > FLT_MIN) {
      spectrum[k] *= filter[k];//Half hann window tappering in favor of high frequencies
      // if(k < N) {
      //   spectrum[N-k] *= filter[k];//Half hann window tappering in favor of high frequencies
      // }
    }
  }
}

//Spectral smoothing with rectangular boxcar or unweighted sliding-average smooth
// This is form Audacity code (HALF SPECTRUM VERSION)
void spectral_smoothing_boxcar(float* spectrum, int kernel_width,int N, bool normalized){
  int k;
  float smoothing_tmp[N+1];
  float t_spectrum[N+1];

  //Initialize smothingbins_tmp
  for (k = 0; k <= N; k++) {
    if (normalized){
      t_spectrum[k] = logf(spectrum[k]);
    }else{
      t_spectrum[k] = spectrum[k];
    }
    //Initialize temporal spectrum
    smoothing_tmp[k] = 0.f;
  }

  for (k = 0; k < N; k++) {
    const int j0 = MAX(0, k - kernel_width);
    const int j1 = MIN(N, k + kernel_width);
    for(int l = j0; l <= j1; ++l) {
      smoothing_tmp[k] += t_spectrum[l];
    }
    smoothing_tmp[k] /= (j1 - j0 + 1);
  }

  for (k = 0; k <= N; k++){
    if (normalized){
      spectrum[k] = expf(smoothing_tmp[k]);
    }else{
      spectrum[k] = smoothing_tmp[k];
    }
  }
}
