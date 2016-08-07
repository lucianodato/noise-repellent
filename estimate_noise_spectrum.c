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

void estimate_noise_spectrum(float* p2,
                             int type_noise_estimation,
                             int fft_size_2,
                             float* noise_spectrum,
                             float* noise_print_min,
      											 float* noise_print_max,
                             int noise_stat_choise){
  int k;
  float overall_max, max_value;

  switch (type_noise_estimation){
    case 1:
      //Manual Capture

      //get min max of the power spectrum
      for(k = 0 ; k <= fft_size_2 ; k++) {
        noise_print_min[k] = MIN(noise_print_min[k], p2[k]);
        noise_print_max[k] = MAX(noise_print_max[k], p2[k]);
      }

      //NOISE SPECTUM COMSTRUCTIOM BASED ON STATISTICS SELECTED

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
      overall_max = max_spectral_value(noise_print_max,fft_size_2);
      max_value = max_spectral_value(noise_spectrum,fft_size_2);

      //Normalize noise spectrum (to get the same level while changing statistic)
      for(k = 0 ; k <= fft_size_2 ; k++) {
        noise_spectrum[k] /= max_value;
        //Rescale it based on overall max_value
        noise_spectrum[k] *= overall_max;
      }

      break;
    case 2:
      //Adaptive noise estimation
      //estimate_noise_loizou(p2,thresh,prev_noise,prev_p2,prev_p_min,prev_speech_p_p)




      break;
  }
}
