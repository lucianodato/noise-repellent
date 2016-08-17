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

void denoise_gain_ps(float over_reduc,
                      int fft_size_2,
                      float* p2,
                      float* noise_thresholds,
                      float* Gk,
                      float* Gk_prev) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (noise_thresholds[k] > FLT_MIN){
      if(p2[k] > FLT_MIN){
        gain = MAX(p2[k]-noise_thresholds[k], 0.f) / p2[k];
      } else {
        gain = 0.f;
      }
      //Use oversustraction
      Fk = over_reduc*(1.0-gain) ;

      if(Fk < 0.0) Fk = 0.0 ;
      if(Fk > 1.0) Fk = 1.0 ;

      Gk_prev[k] = Gk[k];
      Gk[k] =  1.0 - Fk ;

    } else {
      //Otherwise we keep everything as is
      Gk_prev[k] = Gk[k];
      Gk[k] = 1.f;
    }
  } //for
}
void denoise_gain_w(float over_reduc,
                      int fft_size_2,
                      float* p2,
                      float* noise_thresholds,
                      float* Gk,
                      float* Gk_prev) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (noise_thresholds[k] > FLT_MIN){
      float aux = (p2[k] - noise_thresholds[k]);
      if(p2[k] > noise_thresholds[k]){
        gain = aux/(aux+noise_thresholds[k]);
      } else {
        gain = 0.f;
      }
      //Use oversustraction
      Fk = over_reduc*(1.0-gain) ;

      if(Fk < 0.0) Fk = 0.0 ;
      if(Fk > 1.0) Fk = 1.0 ;

      Gk_prev[k] = Gk[k];
      Gk[k] =  1.0 - Fk ;
    } else {
      //Otherwise we keep everything as is
      Gk_prev[k] = Gk[k];
      Gk[k] = 1.f;
    }
  } //for
}

void denoise_gain_mmse(float over_reduc,
                        int option,
                        float alpha_set,
                        int* prev_frame,
                        float* p2,
                        float* p2_prev,
                        float* gain_prev,
                        int fft_size_2,
                        float* Gk,
                        float* Gk_prev,
                        float* noise_thresholds) {
    int k;
    float gain, Fk, Rpost, Rprio, alpha;

    for (k = 0; k <= fft_size_2 ; k++) {
      if (noise_thresholds[k] > FLT_MIN){
        // EM using Wolfe and Godsill optimization
        Rpost = MAX(p2[k]/noise_thresholds[k]-1.f, 0.f);

        if(option == 0){ //Ephraim-Malah
            alpha = alpha_set;
        } else {
          if (Rpost > 0.f) { // CMSR
            alpha = alpha_set; //EM like when Posteriori estimation is null
          } else {
            alpha = 0.f; //Wiener like punctual supression when Posteriori estimation is not null
          }
        }

        if(*(prev_frame) == 1) {
          Rprio = (1.f-alpha)*Rpost + alpha*gain_prev[k]*gain_prev[k]*(p2_prev[k]/noise_thresholds[k]);
        }else{
          Rprio = Rpost;
        }

        //Ephraim-Malah noise suppression, from Godsill and Wolfe 2001 paper (cheaper)
        float r = MAX(Rprio/(1.f+Rprio),FLT_MIN);
        float V = (Rprio/(1.f+Rprio))*(Rpost+1.f);
        gain = sqrtf( r * (1.f+V)/(Rpost+1.f) );

        //Use oversustraction
        Fk = over_reduc*(1.0-gain) ;

        if(Fk < 0.0) Fk = 0.0 ;
        if(Fk > 1.0) Fk = 1.0 ;

        gain =  1.0 - Fk ;

        p2_prev[k] = p2[k];
        gain_prev[k] = gain;
        *(prev_frame) = 1;

        Gk_prev[k] = Gk[k];
        Gk[k] =  gain;

      } else {
        //Otherwise we keep everything as is
        Gk_prev[k] = Gk[k];
        Gk[k] = 1.f;
      }
    } //for
  }

  // void denoise_gain_gss(float alpha,
  //                       float beta,
  //                       float gamma1,
  //                       float gamma2,
  //                       int fft_size_2,
  //                       float* p2,
  //                       float* noise_thresholds,
  //                       float* Gk) {
  //   int k;
  //   float gain, SNRp;
  //   float comp = powf((alpha+beta),2.f);
  //
  //   for (k = 0; k <= fft_size_2 ; k++) {
  //     if (noise_thresholds[k] > FLT_MIN){
  //       SNRp = p2[k]/noise_thresholds[k];
  //       if(SNRp > comp){
  //         gain = powf( (1.f - (alpha / powf(SNRp, (gamma1/2.f) ) ) ) , (gamma2/2.f));
  //       } else {
  //         gain = powf( (beta / powf(SNRp, (gamma1/2.f) ) ) , (gamma2/2.f));
  //       }
  //
  //       Gk[k] = gain;
  //
  //     } else {
  //       //Otherwise we keep everything as is
  //       Gk[k] = 1.f;
  //     }
  //   } //for
  // }
