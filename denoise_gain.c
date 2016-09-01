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


//Spectral Subtraction
void denoise_gain_ss(float over_sustraction,
                      int fft_size_2,
                      float* spectrum,
                      float* noise_thresholds,
                      float* Gk) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (noise_thresholds[k] > FLT_MIN){
      if(spectrum[k] > FLT_MIN){
        gain = MAX(spectrum[k]-noise_thresholds[k], 0.f) / spectrum[k];
      } else {
        gain = 0.f;
      }
      //Use oversustraction
      Fk = over_sustraction*(1.0-gain) ;

      if(Fk < 0.0) Fk = 0.0 ;
      if(Fk > 1.0) Fk = 1.0 ;

      Gk[k] =  1.0 - Fk ;

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for
}
