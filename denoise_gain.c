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


/*Generalized Spectral Subtraction
gamma defines what type of spectral Subtraction is used
  gamma1=gamma2=1 is magnitude substaction
  gamma1=2 gamma2=0.5 is power Subtraction
  gamma1=2 gamma2=1 is wiener filtering
alpha is the oversustraction factor
beta is the spectral flooring factor
reduction_strenght is the other oversustraction designed by the user
so there are 2 oversustraction factors
*/
void denoise_gain_gss(float reduction_strenght,
                      int fft_size_2,
                      float gamma1,
                      float gamma2,
                      float alpha,
                      float beta,
                      float* spectrum,
                      float* noise_thresholds,
                      float* Gk) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (spectrum[k] > FLT_MIN){
      if(powf((noise_thresholds[k]/spectrum[k]),gamma1) < (1.f/(alpha+beta))){
        gain = MAX(powf(1.f-alpha*powf((noise_thresholds[k]/spectrum[k]),gamma1),gamma2),0.f);
      } else {
        gain = MAX(powf(beta*powf((noise_thresholds[k]/spectrum[k]),gamma1),gamma2),0.f);
      }
      //Use reduction_strenght
      Fk = reduction_strenght*(1.f-gain);

      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      Gk[k] =  1.f - Fk;

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for
}

//This version recieves a vector of alphas and betas
void denoise_gain_gss_v(float reduction_strenght,
                      int fft_size_2,
                      float gamma1,
                      float gamma2,
                      float* alpha,
                      float* beta,
                      float* spectrum,
                      float* noise_thresholds,
                      float* Gk) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (spectrum[k] > FLT_MIN){
      if(powf((noise_thresholds[k]/spectrum[k]),gamma1) < (1.f/(alpha[k]+beta[k]))){
        gain = MAX(powf(1.f-alpha[k]*powf((noise_thresholds[k]/spectrum[k]),gamma1),gamma2),0.f);
      } else {
        gain = MAX(powf(beta[k]*powf((noise_thresholds[k]/spectrum[k]),gamma1),gamma2),0.f);
      }
      //Use reduction_strenght
      Fk = reduction_strenght*(1.f-gain);

      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      Gk[k] =  1.f - Fk;

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for
}
