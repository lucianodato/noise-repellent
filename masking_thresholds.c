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

void compute_bark_z(float* bark_z,int fft_size_2, int srate) {
  int k;
  float freq;
  /* compute the bark z value for this frequency bin */
  for(k = 1 ; k <= fft_size_2 ; k++) {
    freq = (float)srate / 2.f /(float)(fft_size_2)*(float)k ; //bin to freq
    bark_z[k] = 7.f*logf(freq/650.f + sqrtf(1.f + (freq/650.f)*(freq/650.f))) ;
  }
}

void compute_johnston_gain(float* bark_z,float** jg_upper,float** jg_lower,int fft_size_2, float tonality_factor) {
  int k,j;
  float bark_diff,johnston,johnston_masked,gain_j;

  //Critical band Analysis - loops take into account bark scale
  for (k = 1; k <= fft_size_2 ; ++k) {
    for(j = k-1 ; j > 0 ; j--) {
      bark_diff = bark_z[k] - bark_z[j];

      //Spreading function
      johnston = 15.81 + 7.5*(bark_diff+0.474) - 17.5*sqrtf(1.f+(bark_diff+0.474)*(bark_diff+0.474));
      //Masking offset
      johnston_masked = johnston - (tonality_factor*(14.5+bark_z[j])+5.5*(1.f - tonality_factor));
      //spread Masking threshold
      gain_j = powf(10.f, johnston_masked/10.f);

      jg_lower[k][k-j] = gain_j;

      if(k - j > 10) break;
    }
    for(j = k ; j <= fft_size_2 ; j++) {
      bark_diff = bark_z[j] - bark_z[k];

      //Spreading function
      johnston = 15.81 + 7.5*(bark_diff+0.474) - 17.5*sqrtf(1.f+(bark_diff+0.474)*(bark_diff+0.474));
      //Masking offset
      johnston_masked = johnston - (tonality_factor*(14.5+bark_z[j])+5.5*(1.f - tonality_factor));
      //spread Masking threshold
      gain_j = powf(10.f, johnston_masked/10.f);

      jg_upper[k][j-k] = gain_j;

      if(j - k > 10) break;
    }
  }
}

void compute_masked(float* spectrum, float* noise_thresholds,float** jg_upper,float** jg_lower, float* masked, int fft_size_2) {
  int j,k;
  float gain_m;

  //renormalization and comparition with absolute threshold
  //Relating the spread masking threshold to the critical band masking thresholds
  for (k = 1; k <= fft_size_2 ; k++) {
    masked[k] = 0.f;

    for(j = k-1 ; j > 0 ; j--) {
      gain_m = jg_lower[k][k-j];

      //take into account absolute threshold of hearing
      if(gain_m < 1.e-2) break;
      if(k - j > 10) break;

      masked[k] += MAX((spectrum[j]-noise_thresholds[j]),0.f)*gain_m;
    }

    for(j = k ; j <= fft_size_2 ; j++) {
      gain_m = jg_upper[k][j-k];

      //take into account absolute threshold of hearing
      if(gain_m < 1.e-2) break;
      if(j - k > 10) break;

      masked[k] += MAX((spectrum[j]-noise_thresholds[j]),0.f)*gain_m;
    }

  }
}
