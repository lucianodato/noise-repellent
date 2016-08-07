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

#include <float.h>
#include <math.h>
#include "estimate_noise_spectrum.c"

static float max_spectral_value(float* noise_print, int N){
  int k;
  float max = 0.f;
  for(k = 0; k <= N; k++){
    if (noise_print[k] > max) max = noise_print[k];
  }
  return max;
}

static float min_spectral_value(float* noise_print, int N){
  int k;
  float min = FLT_MAX;
  for(k = 0; k <= N; k++){
    if (noise_print[k] < min) min = noise_print[k];
  }
  return min;
}

// //Spectral smoothing (Based on Audacity code)
// static void gain_spectral_smoothing(float* gain_spectrum, float* gains, int smoothing_bins,int N){
//   int k;
//   float smoothing_tmp[N+1];
//   int middle_bin = N/4;
//
//   //Initialize smothingbins_tmp
//   for (k = 0; k <= N; ++k) {
//     gains[k] = log(gain_spectrum);
//     smoothing_tmp[k] = 0.f;
//   }
//
//   //do not smooth up to the middle bin not even DC
//   for (k = 0; k < middle_bin; ++k) {
//     smoothing_tmp[k] = gain_spectrum[k];
//   }
//
//   for (k = middle_bin; k < N; ++k) {
//     const int j0 = MAX(middle_bin, k - smoothing_bins);
//     const int j1 = MIN(N, k + smoothing_bins);
//     for(int l = j0; l <= j1; ++l) {
//        smoothing_tmp[k] += gain_spectrum[l];
//     }
//     smoothing_tmp[k] /= (j1 - j0 + 1);
//   }
//
//   for (k = 0; k <= N; ++k){
//     //if (gain_spectrum[k] < 1.f) {
//       gain_spectrum[k] = smoothing_tmp[k];
//   		if(k < N)
//   			gain_spectrum[N-k] = smoothing_tmp[k];
//     //}
//   }
// }


static float gain_weiner(float Yk2, float Dk2) {
  float gain;
  float Xk2 = Yk2 - Dk2;

  if(Yk2 > Dk2)
  gain = (Xk2) / (Xk2+Dk2);
  else
  gain = 0.f;

  return gain;
}

static float gain_power_subtraction(float Yk2, float Dk2) {
  float level = MAX(Yk2-Dk2, 0.f);

  if(Yk2 > FLT_MIN)
  return level/Yk2;
  else
  return 0.f;
}

static float gain_em(float Rprio, float Rpost) {
  float gain;

  //Ephraim-Malah noise suppression, from Godsill and Wolfe 2001 paper (cheaper)
  float r = MAX(Rprio/(1.f+Rprio),FLT_MIN) ;
  float V = (Rprio/(1.f+Rprio))*(Rpost+1.f) ;
  gain = sqrt( r * (1.f+V)/(Rpost+1.f) ) ;

  return gain;
}

void denoise_gain(int denoise_method,
                  float over_reduc,
                  float* p2,
                  float* p2_prev,
                  int fft_size_2,
                  int noise_stat_choise,
                  float* noise_print_min,
                  float* noise_print_max,
                  float* Gk,
                  float* Gk_prev,
                  float* gain_prev,
                  float* noise_spectrum,
                  float alpha_set,
                  int* prev_frame) {
  int k;
  float gain, Fk;

  //----------------------PREPROCESSING-----------------------

  //NOISE SPECTUM COMSTRUCTIOM BASED ON STATISTICS

  //time smoothing for each bin of the captured spectrum
  for(k = 0 ; k <= fft_size_2 ; k++) {
    switch(noise_stat_choise){
      case 0:
      noise_spectrum[k] = noise_print_max[k]; // max spectrum
      break;
      case 1:
      noise_spectrum[k] = noise_print_min[k] + 0.5*(noise_print_max[k] - noise_print_min[k]); // geometric mean spectrum
      break;
    }
  }

  //Find max value of all noise vectors
  float overall_max = max_spectral_value(noise_print_max,fft_size_2);
  float max_value = max_spectral_value(noise_spectrum,fft_size_2);

  //Normalize noise spectrum
  for(k = 0 ; k <= fft_size_2 ; k++) {
    noise_spectrum[k] /= max_value;
    //Rescale it based on overall max_value
    noise_spectrum[k] *= overall_max;
  }

  //Computing gain for selected algorithm
  for (k = 0; k <= fft_size_2 ; k++) {
    gain = 0.f;
    if (noise_spectrum[k] > FLT_MIN){
      //We can compute gain if print was previously captured
      switch (denoise_method) {// supression rule
        case 0: // Wiener Filter
          gain = gain_weiner(p2[k], noise_spectrum[k]) ;
          break;
        case 1: // Power Subtraction
          gain = gain_power_subtraction(p2[k], noise_spectrum[k]) ;
          break;
        case 2:
          // Ephraim-Mallat - Using CMSR rule
          float Rpost = MAX(p2[k]/noise_spectrum[k]-1.f, 0.f);

          float alpha;
          if (Rpost > 0.f){ // Canazza-Mian Condition
            alpha = alpha_set; // Traditional EM when Rpost is high
          }else{
            alpha = 0.f; //Wiener like when low Rpost
          }

          float Rprio;

          if(*(prev_frame) == 1) {
            Rprio = (1.f-alpha)*Rpost + alpha*gain_prev[k]*gain_prev[k]*(p2_prev[k]/noise_spectrum[k]);
          }else{
            Rprio = Rpost;
          }

          gain = gain_em(Rprio, Rpost);

          p2_prev[k] = p2[k];
          gain_prev[k] = gain;
          *(prev_frame) = 1;
          break;
      }

      //To avoid excesive distortion limit gain to be applied
      Fk = over_reduc*(1.f-gain);

      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      Gk[k] =  1.f - Fk;
      Gk_prev[k] = Gk[k];

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for

}
