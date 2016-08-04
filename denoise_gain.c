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
                  float amount,
                  float* p2,
                  float* p2_prev,
                  int fft_size_2,
                  float* Gk,
                  float* gain_prev,
                  float* noise_spectrum,
                  float alpha_set,
                  int* prev_frame) {
  int k;
  float gain, Fk;

  //Precalculation for Canazza-Mian Rule
  float rpost_sum = 0.f;
  if (denoise_method == 2){
    for (k = 0; k <= fft_size_2 ; k++) {
      rpost_sum += MAX(p2[k]/noise_spectrum[k]-1.f, 0.f);
    }
  }

  //Computing gain and applying the Reduction
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
          if (Rpost > 0.f && rpost_sum < 0.f){ // Canazza-Mian Condition
            alpha = 0.f; // Wiener like
          }else{
            alpha = alpha_set; // Traditional EM
          }
          float Rprio;

          if(*(prev_frame) == 1) {
            Rprio = (1.f-alpha)*Rpost+alpha*gain_prev[k]*gain_prev[k]*p2_prev[k]/noise_spectrum[k];
          }else{
            Rprio = Rpost;
          }

          gain = gain_em(Rprio, Rpost);

          p2_prev[k] = p2[k];
          gain_prev[k] = gain;
          *(prev_frame) = 1;
          break;
      }

      Fk = amount*(1.f-gain);

      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      Gk[k] =  1.f - Fk;
    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for
}
