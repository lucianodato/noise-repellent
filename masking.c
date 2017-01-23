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
#define ALPHA_MAX 6.f
#define ALPHA_MIN 1.f
#define BETA_MAX 0.02f
#define BETA_MIN 0.0f
#define BINS_PER_BARK_BAND_LIMIT 10
#define NOISE_FLOOR 1.e-2f


//fft to bark bilinear transform
void compute_bark_z(float* bark_z,int fft_size_2, int srate) {
  int k;
  float freq;
  /* compute the bark z value for this frequency bin */
  for(k = 1 ; k <= fft_size_2 ; k++) {
    freq = (float)srate / 2.f /(float)(fft_size_2)*(float)k ; //bin to freq
    bark_z[k] = 7.f*logf(freq/650.f + sqrtf(1.f + (freq/650.f)*(freq/650.f))) ;
  }
}

void compute_masking_thresholds(float* bark_z,
                                float* spectrum,
                                int fft_size_2,
                                float* masking_thresholds,
                                float tonality_factor) {
  int k,j;
  float bark_diff,johnston,johnston_masked,gain_j;

  //Critical band Analysis
  for (k = 1; k <= fft_size_2 ; ++k) {

    //Lower frequencies bellow k
    for(j = k-1 ; j > 0 ; j--) {
      //Limit number of bins to save processing
      if(k - j > BINS_PER_BARK_BAND_LIMIT) break;

      //frequency range
      bark_diff = bark_z[k] - bark_z[j];

      //Spreading function
      johnston = 15.81 + 7.5*(bark_diff+0.474) - 17.5*sqrtf(1.f+(bark_diff+0.474)*(bark_diff+0.474));
      //Masking offset
      johnston_masked = johnston - (tonality_factor*(14.5+bark_z[j])+5.5*(1.f - tonality_factor));
      //spread Masking threshold
      gain_j = powf(10.f, johnston_masked/10.f);

      //take into account absolute threshold of hearing
      if(gain_j < NOISE_FLOOR) break;

      //Relating the spread masking threshold to the critical band masking thresholds
      masking_thresholds[k] += spectrum[j]*gain_j;
    }

    //Upper frequencies above k
    for(j = k ; j <= fft_size_2 ; j++) {
      //Limit number of bins to save processing
      if(j - k > BINS_PER_BARK_BAND_LIMIT) break;

      //frequency range
      bark_diff = bark_z[j] - bark_z[k];

      //Spreading function
      johnston = 15.81 + 7.5*(bark_diff+0.474) - 17.5*sqrtf(1.f+(bark_diff+0.474)*(bark_diff+0.474));
      //Masking offset
      johnston_masked = johnston - (tonality_factor*(14.5+bark_z[j])+5.5*(1.f - tonality_factor));
      //spread Masking threshold
      gain_j = powf(10.f, johnston_masked/10.f);

      //take into account absolute threshold of hearing
      if(gain_j < NOISE_FLOOR) break;

      //Relating the spread masking threshold to the critical band masking thresholds
      masking_thresholds[k] += spectrum[j]*gain_j;
    }
  }
}

//alpha and beta computation to be used in general spectral Sustraction
void compute_alpha_and_beta(float* bark_z,
                            float* fft_p2,
                            float* noise_thresholds_p2,
                            int fft_size_2,
                            float* alpha,
                            //float* beta,
                            float* max_masked,
                            float* min_masked,
                            float masking) {

  int k;
  float masking_thresholds[fft_size_2+1];
  float estimated_clean[fft_size_2+1];
  float SFM,tonality_factor,weight1,sum_p = 0.f,sum_log_p = 0.f;	//,weight2;

  //Noise masking threshold must be computed from a clean signal
  //therefor we aproximate a clean signal using a power Sustraction over
  //the original noisy one

  //basic spectral Sustraction to estimate clean signal
  for (k = 0; k <= fft_size_2; k++) {
    estimated_clean[k] =  MAX(fft_p2[k]-noise_thresholds_p2[k],0.f);
    //For Geometric and Arithmetic mean computing
    sum_p += estimated_clean[k];
    sum_log_p += log10f(estimated_clean[k]);
  }

  //spectral flatness measure using Geometric and Arithmetic means of the spectrum cleaned previously
  SFM = 10.f*((sum_log_p/(float)fft_size_2) - log10f(sum_p/(float)fft_size_2));//this value is in db scale

  //Tonality factor in db scale
  tonality_factor = MIN(SFM/-60.f, 1.f);


  /*Now we can compute noise masking threshold from this clean signal
   *1- we need to compute the energy in the bark scale
   *2- Convolution with a spreading function
   *3- Sustraction of the offset depending of noise masking tone masking
   *4- renormalization and comparition with the absolute threshold of hearing
   */
  compute_masking_thresholds(bark_z,
                             estimated_clean,
                             fft_size_2,
                             masking_thresholds,
                             tonality_factor);

  /*Get alpha and beta based on masking thresholds
  *beta and alpha values would adapt based on masking thresholds
  *frame to frame for optimal oversustraction and noise floor parameter in each one
  *noise floor would better be controled by user using the amount of reduction
  *so beta is not modified
  */

  //First we need the maximun and the minimun value of the masking threshold
  *(max_masked) = MAX(max_spectral_value(masking_thresholds,fft_size_2),*(max_masked));
  *(min_masked) = MIN(min_spectral_value(masking_thresholds,fft_size_2),*(min_masked));


  weight1 = (masking - ALPHA_MIN)/(*(min_masked) - *(max_masked));
  //float weight2 = (BETA_MAX - BETA_MIN)/(*(min_masked) - *(max_masked)));

  for (k = 0; k <= fft_size_2; k++) {
    //new alpha and beta vector
    if(masking_thresholds[k] == *(max_masked)){
       alpha[k] = ALPHA_MIN;
       //beta[k] = BETA_MIN;
    }
    if(masking_thresholds[k] == *(min_masked)){
       alpha[k] = masking;
       //beta[k] = BETA_MAX;
    }
    if(masking_thresholds[k] < *(max_masked) && masking_thresholds[k] > *(min_masked)){
       //Linear interpolation of the value between max and min masked threshold values
       alpha[k] = ALPHA_MIN + weight1 * (masking_thresholds[k]- ALPHA_MIN);
       //beta[k] = BETA_MIN + weight2 * (masking_thresholds[k]- BETA_MIN);
    }
  }
}
