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

//masking thresholds values recomended by virag
#define ALPHA_MAX 10.0
#define ALPHA_MIN 1.0
#define BETA_MAX 0.02
#define BETA_MIN 0.0

void compute_bark_z(float* bark_z,int fft_size_2, int srate) {
  int k;
  float freq;
  /* compute the bark z value for this frequency bin */
  for(k = 1 ; k <= fft_size_2 ; k++) {
    freq = (float)srate / 2.f /(float)(fft_size_2)*(float)k ; //bin to freq
    bark_z[k] = 7.f*logf(freq/650.f + sqrtf(1.f + (freq/650.f)*(freq/650.f))) ;
  }
}

void compute_masking_thresholds(float* bark_z,float* spectrum, float* noise_thresholds,int fft_size_2,float* masked, float tonality_factor) {
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

      //take into account absolute threshold of hearing
      if(gain_j < 1.e-2) break;
      if(k - j > 10) break;

      //Relating the spread masking threshold to the critical band masking thresholds
      masked[k] += spectrum[j]*gain_j;
    }
    for(j = k ; j <= fft_size_2 ; j++) {
      bark_diff = bark_z[j] - bark_z[k];

      //Spreading function
      johnston = 15.81 + 7.5*(bark_diff+0.474) - 17.5*sqrtf(1.f+(bark_diff+0.474)*(bark_diff+0.474));
      //Masking offset
      johnston_masked = johnston - (tonality_factor*(14.5+bark_z[j])+5.5*(1.f - tonality_factor));
      //spread Masking threshold
      gain_j = powf(10.f, johnston_masked/10.f);

      //take into account absolute threshold of hearing
      if(gain_j < 1.e-2) break;
      if(j - k > 10) break;

      //Relating the spread masking threshold to the critical band masking thresholds
      masked[k] += spectrum[j]*gain_j;
    }
  }
}

//alpha and beta computation to be used in general spectral Sustraction
void compute_alpha_and_beta(float* bark_z,
                            float* fft_p2,
                            float* fft_magnitude,
                            float* noise_thresholds,
                            int fft_size_2,
                            float* alpha,
                            //float* beta,
                            float max_masked,
                            float min_masked) {

  int k;
  float estimated_clean[fft_size_2+1];
  float masked[fft_size_2+1];
  float gmean_value,mean_value,SFM,tonality_factor,weight1;//,weight2

  //Noise masking threshold must be computed from a clean signal
  //therefor we aproximate a clean signal using a power Sustraction over
  //the original noisy one

  //basic power Sustraction to estimate clean signal
  for (k = 0; k <= fft_size_2; k++) {
    estimated_clean[k] =  MAX(fft_p2[k]-powf(noise_thresholds[k],2),0.f);
  }

  //spectral flatness measure using Geometric and Arithmetic means of the spectrum cleaned previously
  gmean_value = spectral_gmean(fft_size_2,estimated_clean);
  mean_value = spectral_mean(fft_size_2,estimated_clean);
  SFM = 10.f*log10f(gmean_value/mean_value);//this value is in db scale

  //Tonality factor in db scale
  tonality_factor = MIN(SFM/60.f, 1.f);

  /*Now we can compute noise masking threshold from this clean signal
   *1- we need to compute the energy in the bark scale
   *2- Convolution with a spreading function
   *3- Sustraction of the offset depending of noise masking tone masking
   *4- renormalization and comparition with the absolute threshold of hearing
   */
  compute_masking_thresholds(bark_z,
                             estimated_clean,
                             noise_thresholds,
                             fft_size_2,
                             masked,
                             tonality_factor);

  /*Get alpha and beta based on masking thresholds
   *beta and alpha values would adapt based on masking thresholds
   *frame to frame for optimal oversustraction and noise floor parameter in each one
   *noise floor would better be controled by user using the amount of reduction
   *so beta is not modified
  */

  //First we need the maximun and the minimun value of the masking threshold
  max_masked = MAX(max_spectral_value(masked,fft_size_2),max_masked);
  min_masked = MIN(min_spectral_value(masked,fft_size_2),min_masked);
  weight1 = (ALPHA_MAX - ALPHA_MIN)/(min_masked - max_masked);
  //float weight2 = (BETA_MAX - BETA_MIN)/(min_masked - max_masked);

  for (k = 0; k <= fft_size_2; k++) {
    //new alpha and beta vector
    if(masked[k] == max_masked){
        alpha[k] = ALPHA_MIN;
        //beta[k] = BETA_MIN;
    }
    if(masked[k] == min_masked){
        alpha[k] = ALPHA_MAX;
        //beta[k] = BETA_MAX;
    }
    if(masked[k] < max_masked && masked[k] > min_masked){
        //Linear interpolation of the value between max and min masked threshold values
        alpha[k] = ALPHA_MIN + weight1 * (masked[k]- ALPHA_MIN);
        //beta[k] = BETA_MIN + weight2 * (masked[k]- BETA_MIN);
    }
  }
}
