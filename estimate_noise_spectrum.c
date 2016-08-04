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

void estimate_noise_spectrum(int noise_stat_choise,
                             float* p2,
                             int type_noise_estimation,
                             float* noise_print_min,
                             float* noise_print_max,
                             float* noise_print_avg,
                             int count,
                             int fft_size_2,
                             float* noise_spectrum){
  int k;
  switch (type_noise_estimation){
    case 1:
      //Manual Capture

      //Increment number of windows processed to correctly average
      count++;

      //get min max and average of thes power spectrum
      for(k = 0 ; k <= fft_size_2 ; k++) {
        noise_print_min[k] = MIN(noise_print_min[k], p2[k]);
        noise_print_max[k] = MAX(noise_print_max[k], p2[k]);
        noise_print_avg[k] += (p2[k]-noise_print_avg[k])/count; //Moving average
      }

      //time smoothing for each bin of the captured spectrum
      for(k = 0 ; k <= fft_size_2 ; k++) {
        switch(noise_stat_choise){
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
      }
      break;
    case 2:
      //Auto Capture
      break;
  }
}
