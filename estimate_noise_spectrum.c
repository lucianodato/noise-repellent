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

#define MANUAL 1
#define ADAPTIVE 2

//For louizou algorithm
#define N_SMOOTH 0.7 //Smoothing over the power spectrum [1 - previous / 0 - actual]
#define BETA 0.8 //Adaption time of the local minimun [1 - slower / 0 - faster]
#define GAMMA 0.998 //Smoothing factor over local minimun [1 - previous / 0 - actual]
#define ALPHA_P 0.2 //smoothing constant over speech presence [1 - previous / 0 - actual]
#define ALPHA_D 0.85 //time–frequency dependent smoothing [0-1] [1 - previous / 0 - actual]

//For exponential smoothing of the time smooth noise print
#define TVEL 0.1

static void estimate_noise_loizou(float thresh,
                      int fft_size_2,
                      float* p2,
                      float* s_pow_spec,
                      float* prev_s_pow_spec,
                      float* a_noise_spectrum,
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
    //1- Compute the noisy speech power spectrum
    s_pow_spec[k] = N_SMOOTH * prev_s_pow_spec[k] + (1.f-N_SMOOTH) * p2[k]; //interpolation between

    //2- Compute the local minimum of noisy speech
    if(prev_p_min < s_pow_spec) {
      p_min[k] = GAMMA * prev_p_min[k] + ((1.f-GAMMA)/(1.f-BETA)) * (s_pow_spec[k] - BETA * prev_s_pow_spec[k]);
    } else {
      p_min[k] = s_pow_spec[k];
    }

    //3- Compute ratio of noisy speech power spectrum to its local minimum
    ratio_ns = s_pow_spec[k]/p_min[k];

    //4- Compute the indicator function I for speech present/absent detection
    if(ratio_ns > thresh) { //thresh could be freq dependant
      speech_p_d[k] = 1.f; //present
    } else {
      speech_p_d[k] = 0.f; //absent
    }

    //5- Calculate speech presence probability using first-order recursion
    speech_p_p[k] = ALPHA_P * prev_speech_p_p[k] + (1.f-ALPHA_P) * speech_p_d[k];

    //6- Compute time-frequency dependent smoothing constant
    freq_s[k] = ALPHA_D + (1.f-ALPHA_D) * speech_p_p[k];

    //7- Update noise estimate D using time-frequency dependent smoothing factor α s (λ,k).
    a_noise_spectrum[k] = freq_s[k] * prev_noise[k] + (1.f-freq_s[k]) * p2[k];
  }
}

//Adaptive noise estimation
void auto_capture_noise(float* p2,
                        int fft_size_2,
                        float* a_noise_spectrum,
                        float thresh,
                        float* prev_a_noise,
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
                        a_noise_spectrum,
                        prev_a_noise,
                        p_min,
                        prev_p_min,
                        speech_p_p,
                        prev_speech_p_p);

  //Update previous variables
  for(k = 0 ; k <= fft_size_2 ; k++) {
    prev_a_noise[k] = a_noise_spectrum[k];
    prev_s_pow_spec[k] = s_pow_spec[k];
    prev_p_min[k] = p_min[k];
    prev_speech_p_p[k] = speech_p_p[k];
  }
}

void get_noise_statistics(float* p2,
                         int fft_size_2,
                         float* noise_print_min,
  											 float* noise_print_max,
                         float* noise_print_g_mean,
                         float* noise_print_average,
                         float* n_window_count) {
  int k;

  *(n_window_count) += 1.f;

  //Get noise time statistics
  for(k = 0 ; k <= fft_size_2 ; k++) {
    noise_print_min[k] = MIN(noise_print_min[k], p2[k]);
    noise_print_max[k] = MAX(noise_print_max[k], p2[k]);
    if(*(n_window_count) == 1){
      noise_print_average[k] = p2[k];
    } else {
      noise_print_average[k] += ((p2[k] - noise_print_average[k])/ *(n_window_count));
    }
  }

  //Finish geometric mean calculations
  for(k = 0 ; k <= fft_size_2 ; k++) {
    noise_print_g_mean[k] = noise_print_min[k] + 0.5*(noise_print_max[k] - noise_print_min[k]);
  }
}

//Manual Capture
void estimate_noise_thresholds(int fft_size_2,
                             int stat_choise,
                             float* noise_thresholds,
                             float* noise_print_max,
      											 float* noise_print_g_mean,
                             float* noise_print_average){
  int k;
  //NOISE THRESHOLDS SPECTUM CONSTRUCTIOM BASED ON STATISTIC SELECTED

  //Do some more stuff to decide what is the threshold of each band
  //Like amplifying the definitive spectrum to catch more noise_spectrum
  //sensitivity control etc

  //Temporary
  for(k = 0 ; k <= fft_size_2 ; k++) {
    switch (stat_choise) {
      case 0:
      noise_thresholds[k] = noise_print_max[k];
      break;
      case 1:
      noise_thresholds[k] = noise_print_g_mean[k];
      break;
      case 2:
      noise_thresholds[k] = noise_print_average[k];
      break;
    }
  }

}
