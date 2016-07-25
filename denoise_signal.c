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

 static float gain_weiner(float Yk2, float Dk2)
 {
   float gain;
   float Xk2 = Yk2 - Dk2;

   if(Yk2 > Dk2)
    gain = (Xk2) / (Xk2+Dk2);
   else
    gain = 0.0;

   return gain;
 }

 static float gain_power_subtraction(float Yk2, float Dk2)
 {
   float level = MAX(Yk2-Dk2, 0.0);

   if(Yk2 > FLT_MIN)
    return level/Yk2;
   else
    return 0.0;
 }

static void denoise_signal(int noise_mean_choise,int denoise_method,float amount,float* spectrum,float* noise_print_min,float* noise_print_max,float* noise_print_avg,int fft_size_2) {

  float noise_spectrum[fft_size_2];
  float Y2[fft_size_2];
  int k;

  //convert input to power spectrum
  for(k = 1 ; k <= fft_size_2 ; k++) {
    switch(noise_mean_choise){
      case 0:
        noise_spectrum[k] = noise_print_max[k]; // max spectrum
        break;
      case 1:
        noise_spectrum[k] = noise_print_min[k] + 0.5*(noise_print_max[k] - noise_print_min[k]); // geometric mean spectrum
        break;
      case 2:
        noise_spectrum[k] = noise_print_avg[k]; // mean spectrum
        break;
    }
    Y2[k] = spectrum[k]*spectrum[k]; //Signal power spectrum
  }

  //Computing gain and applying the Reduction
  for (k = 1; k <= fft_size_2 ; k++) {
    if (noise_spectrum[k] > FLT_MIN) {
      float gain = 0, Fk, Gk;

      switch (denoise_method) {// supression rule
        case 0: // Wiener Filter
          gain = gain_weiner(Y2[k], noise_spectrum[k]) ;
          break;
        case 1: // Power Subtraction
          gain = gain_power_subtraction(Y2[k], noise_spectrum[k]) ;
          break;
      }

      Fk = amount*(1.0-gain);

      if(Fk < 0.0) Fk = 0.0;
      if(Fk > 1.0) Fk = 1.0;

      Gk =  1.0 - Fk;

      spectrum[k] = spectrum[k] * Gk;

    }
  }
}
