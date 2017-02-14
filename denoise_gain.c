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

#define SCALING_FACTOR 10.f //scaling factor for non linear power sustraction

//Non linear Power Sustraction
void nonlinear_power_sustraction(float reduction_strenght,
                       float snr_influence,
                       int fft_size_2,
                       float* spectrum,
                       float* noise_thresholds,
                       float* Gk) {
  int k;
  float gain, Fk, alpha;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (noise_thresholds[k] > FLT_MIN){
      if(spectrum[k] > 0.f){
        if(snr_influence > 0){
          alpha = snr_influence + sqrtf(spectrum[k]/noise_thresholds[k]);
        }else{
          alpha = 1.f;//Non linear spectral sustraction off
        }
        gain = MAX(spectrum[k]-alpha*noise_thresholds[k], 0.f) / spectrum[k];
      } else {
        gain = 0.f;
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
  }
}
