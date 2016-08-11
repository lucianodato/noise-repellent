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

void post_filter (float omega,float* p2, float* Gk, float SNR_thresh, int fft_size_2) {
  float SNR_estim, num = 0.f, den = 0.f, N;
  float post_filter[fft_size_2+1];
  int k;

  //Estimate Post-Filter
  for (k = 0; k <= fft_size_2 ; k++) {
    num += powf(Gk[k]*p2[k],2);
    den += powf(p2[k],2);
    if (k<fft_size_2){
      num += powf(Gk[k]*p2[fft_size_2-k],2);
      den += powf(p2[fft_size_2-k],2);
    }
  }

  SNR_estim = num/(den+FLT_MIN);

  if (SNR_estim >= SNR_thresh){
    N = 1.f;
  }else{
    N = 2.f*(1.f+(nearbyintf(SNR_estim/SNR_thresh)*omega)) + 1.f;
  }

  for (k = 0; k <= fft_size_2 ; k++) {
    if(k<N){
      post_filter[k] = 1.f/N;
    }else{
      post_filter[k] = 0.f;
    }
    //Apply Post-Filter
    Gk[k] *= post_filter[k];
  }
}

void denoise_gain(float over_reduc,
                  float* p2,
                  int fft_size_2,
                  float* Gk,
                  float* noise_thresholds) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (noise_thresholds[k] > FLT_MIN){
      // Power Subtraction
      if(p2[k] > FLT_MIN)
        gain = MAX(p2[k]-noise_thresholds[k], 0.f)/p2[k];
      else
        gain = 0.f;

      //Apply over sustraction
      Fk = over_reduc*(1.f-gain);

      //Limit gain to be applied
      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      gain =  1.f - Fk;

      //Assing gain to gain spectrum
      Gk[k] =  gain;

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for
}
