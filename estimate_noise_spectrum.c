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

#include "extra_functions.c"

//For louizou algorithm
#define N_SMOOTH 0.7 //Smoothing over the power spectrum [0.9 - previous / 0.7 - actual]
#define BETA_AT 0.8 //Adaption time of the local minimun [1 - slower / 0 - faster]
#define GAMMA 0.998 //Smoothing factor over local minimun [1 - previous / 0 - actual]
#define ALPHA_P 0.2 //smoothing constant over speech presence [1 - previous / 0 - actual]
#define ALPHA_D 0.99 //time–frequency dependent smoothing [0-1] [1 - previous / 0 - actual]

static void estimate_noise_loizou(float* thresh,
                      int fft_size_2,
                      float* p2,
                      float* s_pow_spec,
                      float* prev_s_pow_spec,
                      float* noise_thresholds,
                      float* prev_noise,
                      float* p_min,
                      float* prev_p_min,
                      float* speech_p_p,
                      float* prev_speech_p_p) {

  int k;
  float ratio_ns = 0.f;
  float freq_s[fft_size_2+1];
  float speech_p_d[fft_size_2+1];

  for(k = 0 ; k <= fft_size_2 ; k++) {
    //1- Smooth between current and past noisy speech power spectrum
    s_pow_spec[k] = N_SMOOTH * prev_s_pow_spec[k] + (1.f-N_SMOOTH) * p2[k]; //interpolation between

    //2- Compute the local minimum of noisy speech
    if(prev_p_min < s_pow_spec) {
      p_min[k] = GAMMA * prev_p_min[k] + ((1.f-GAMMA)/(1.f-BETA_AT)) * (s_pow_spec[k] - BETA_AT * prev_s_pow_spec[k]);
    } else {
      p_min[k] = s_pow_spec[k];
    }

    //3- Compute ratio of noisy speech power spectrum to its local minimum
    ratio_ns = s_pow_spec[k]/p_min[k];

    //4- Compute the indicator function I for speech present/absent detection
    if(ratio_ns > thresh[k]) { //thresh could be freq dependant (it is not a neasure related to dB)
      speech_p_d[k] = 1.f; //present
    } else {
      speech_p_d[k] = 0.f; //absent
    }

    //5- Calculate speech presence probability using first-order recursion
    speech_p_p[k] = ALPHA_P * prev_speech_p_p[k] + (1.f-ALPHA_P) * speech_p_d[k];

    //6- Compute time-frequency dependent smoothing constant
    freq_s[k] = ALPHA_D + (1.f-ALPHA_D) * speech_p_p[k];

    //7- Update noise estimate D using time-frequency dependent smoothing factor α s (λ,k).
    noise_thresholds[k] = freq_s[k] * prev_noise[k] + (1.f-freq_s[k]) * p2[k];
  }
}


//Automatic noise threshold estimation
void auto_capture_noise(float* p2,
                        int fft_size_2,
                        float* noise_thresholds,
                        float* thresh,
                        float* prev_noise_thresholds,
                        float* s_pow_spec,
                        float* prev_s_pow_spec,
                        float* p_min,
                        float* prev_p_min,
                        float* speech_p_p,
                        float* prev_speech_p_p){
  int k;

  //Loizou noise-estimation algorithm for highly non-stationary environments
  estimate_noise_loizou(thresh,
                        fft_size_2,
                        p2,
                        s_pow_spec,
                        prev_s_pow_spec,
                        noise_thresholds,
                        prev_noise_thresholds,
                        p_min,
                        prev_p_min,
                        speech_p_p,
                        prev_speech_p_p);

  //Update previous variables
  for(k = 0 ; k <= fft_size_2 ; k++) {
    prev_noise_thresholds[k] = noise_thresholds[k];
    prev_s_pow_spec[k] = s_pow_spec[k];
    prev_p_min[k] = p_min[k];
    prev_speech_p_p[k] = speech_p_p[k];
  }

  //noise_thresholds should be a magnitude spectrum as
  //general spectral sustraction recieves magnitude spectrum
  for (k = 0; k <= fft_size_2; k++) {
    noise_thresholds[k] = sqrtf(noise_thresholds[k]);
  }
}

//Manual Capture threshold estimation
void get_noise_statistics(float* spectrum,
                         int fft_size_2,
                         float* noise_thresholds,
                         float* window_count) {
  int k;

  *(window_count) += 1.f;

  //Get noise thresholds based on averageing the input noise signal between frames
  for(k = 0 ; k <= fft_size_2 ; k++) {
    if(*(window_count) == 1){
      noise_thresholds[k] = spectrum[k];
    } else {
      noise_thresholds[k] += ((spectrum[k] - noise_thresholds[k])/ *(window_count)); //rolling mean
    }
  }
}
