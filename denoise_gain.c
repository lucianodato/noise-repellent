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

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

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

void denoise_gain(int denoise_method,bool use_mag,float amount,float* mag,float* p2,int fft_size_2,float* Gk,float* noise_spectrum) {
  int k;
  float gain, Fk;

  //Computing gain and applying the Reduction
  for (k = 0; k <= fft_size_2 ; k++) {
    gain = 0;

    if (use_mag){
      switch (denoise_method) {// supression rule
        case 0: // Wiener Filter
        gain = gain_weiner(mag[k], noise_spectrum[k]) ;
        break;
        case 1: // Power Subtraction
        gain = gain_power_subtraction(mag[k], noise_spectrum[k]) ;
        break;
      }
    }else{
      switch (denoise_method) {// supression rule
        case 0: // Wiener Filter
        gain = gain_weiner(p2[k], noise_spectrum[k]) ;
        break;
        case 1: // Power Subtraction
        gain = gain_power_subtraction(p2[k], noise_spectrum[k]) ;
        break;
      }
    }

    Fk = amount*(1.f-gain);

    if(Fk < 0.f) Fk = 0.f;
    if(Fk > 1.f) Fk = 1.f;

    Gk[k] =  1.f - Fk;
  }

}
