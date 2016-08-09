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

void compute_bark_z(float* bark_z,int fft_size_2, int srate) {
  int k;
  float freq;
  /* compute the bark z value for this frequency bin */
  for(k = 1 ; k <= fft_size_2 ; k++) {
    freq = (float)srate / 2.f /(float)(fft_size_2)*(float)k ;
    bark_z[k] = 7.f*logf(freq/650.f + sqrtf(1.f + (freq/650.f)*(freq/650.f))) ;
  }
}

void compute_johnston_gain(float* bark_z,float** jg_upper,float** jg_lower,int fft_size_2, float tonality_factor) {
  int k,j;
  float bark_diff,johnston,johnston_masked,gain_j;

  for (k = 1; k <= fft_size_2 ; ++k) {
    for(j = k-1 ; j > 0 ; j--) {
      bark_diff = bark_z[k] - bark_z[j];

      johnston = 15.81 + 7.5*(bark_diff+0.474) - 17.5*sqrtf(1.f+(bark_diff+0.474)*(bark_diff+0.474));
      johnston_masked = johnston - (tonality_factor*(14.5+bark_z[j])+5.5*(1.f - tonality_factor));
      gain_j = powf(10.f, johnston_masked/10.f);

      jg_lower[k][k-j] = gain_j;

      if(k - j > 10) break;
    }
    for(j = k ; j <= fft_size_2 ; j++) {
      bark_diff = bark_z[j] - bark_z[k];

      johnston = 15.81 + 7.5*(bark_diff+0.474) - 17.5*sqrtf(1.f+(bark_diff+0.474)*(bark_diff+0.474));
      johnston_masked = johnston - (tonality_factor*(14.5+bark_z[j])+5.5*(1.f - tonality_factor));
      gain_j = powf(10.f, johnston_masked/10.f);

      jg_upper[k][j-k] = gain_j;

      if(j - k > 10) break;
    }
  }
}

void compute_masked(float* p2, float* noise_spectrum,float** jg_upper,float** jg_lower, float* masked, int fft_size_2) {
  int j,k;
  float gain_m;

  for (k = 1; k <= fft_size_2 ; k++) {
    masked[k] = 0.f;

    for(j = k-1 ; j > 0 ; j--) {
      gain_m = jg_lower[k][k-j];
      if(k - j > 10) break;

      masked[k] += MAX((p2[j]-noise_spectrum[j]),0.f)*gain_m;
    }

    for(j = k ; j <= fft_size_2 ; j++) {
      gain_m = jg_upper[k][j-k];

      if(gain_m < 1.e-2) break;
      if(j - k > 10) break;

      masked[k] += MAX((p2[j]-noise_spectrum[j]),0.f)*gain_m;
    }
  }
}

static float gain_weiner(float Yk2, float Dk2) {
  float gain_w;
  float Xk2 = Yk2 - Dk2;

  if(Yk2 > Dk2)
  gain_w = (Xk2) / (Xk2+Dk2);
  else
  gain_w = 0.f;

  return gain_w;
}

static float gain_power_subtraction(float Yk2, float Dk2) {
  float level = MAX(Yk2-Dk2, 0.f);

  if(Yk2 > FLT_MIN)
  return level/Yk2;
  else
  return 0.f;
}

static float gain_em(float Rprio, float Rpost) {
  float gain_em;

  //Ephraim-Malah noise suppression, from Godsill and Wolfe 2001 paper (cheaper)
  float r = MAX(Rprio/(1.f+Rprio),FLT_MIN);
  float V = (Rprio/(1.f+Rprio))*(Rpost+1.f);
  gain_em = sqrtf( r * (1.f+V)/(Rpost+1.f) );

  return gain_em;
}

void denoise_gain(float denoise_method,
                  float over_reduc,
                  float* p2,
                  float* p2_prev,
                  int fft_size_2,
                  float* Gk,
                  float* Gk_prev,
                  float* gain_prev,
                  float* noise_spectrum,
                  float alpha_set,
                  int* prev_frame,
                  float* masked,
                  float** jg_upper,
                  float** jg_lower) {
  int k;
  float gain, Fk, Rpost = 0.f, Rprio  = 0.f, alpha  = 0.f;

  if (denoise_method == 4.f) {
    compute_masked(p2,noise_spectrum,jg_upper,jg_lower,masked,fft_size_2);
  }

  //Computing gain for selected algorithm
  for (k = 0; k <= fft_size_2 ; k++) {

    gain = 0.f;
    Fk = 0.f;
    Rpost = 0.f;
    Rprio  = 0.f;
    alpha  = 0.f;

    if (noise_spectrum[k] > FLT_MIN){
      //We can compute gain if print was previously captured
      switch ((int)denoise_method) {// supression rule
        case 0: // Wiener Filter
          gain = gain_weiner(p2[k], noise_spectrum[k]) ;
          break;
        case 1: // Power Subtraction
          gain = gain_power_subtraction(p2[k], noise_spectrum[k]) ;
          break;
        case 2:
          // Ephraim-Mallat
          Rpost = MAX(p2[k]/noise_spectrum[k]-1.f, 0.f);

          if(*(prev_frame) == 1) {
            Rprio = (1.f-alpha_set)*Rpost + alpha_set*gain_prev[k]*gain_prev[k]*(p2_prev[k]/noise_spectrum[k]);
          }else{
            Rprio = Rpost;
          }

          gain = gain_em(Rprio, Rpost);

          p2_prev[k] = p2[k];
          gain_prev[k] = gain;
          *(prev_frame) = 1;
          break;
        case 3:
          // CMSR (modified EM)
          Rpost = MAX(p2[k]/noise_spectrum[k]-1.f, 0.f);

          if (Rpost > 0.f) {
            alpha = alpha_set; //EM like when Posteriori estimation is null
          } else {
            alpha = 0.f; //Wiener like punctual supression when Posteriori estimation is not null
          }

          if(*(prev_frame) == 1) {
            Rprio = (1.f-alpha)*Rpost + alpha*gain_prev[k]*gain_prev[k]*(p2_prev[k]/noise_spectrum[k]);
          }else{
            Rprio = Rpost;
          }

          gain = gain_em(Rprio, Rpost);

          p2_prev[k] = p2[k];
          gain_prev[k] = gain;
          *(prev_frame) = 1;
          break;
        case 4:
          // CMSR (modified using masking thresholds as noise suppression rule)
          Rpost = MAX(p2[k]/noise_spectrum[k]-1.f, 0.f);

          if (Rpost > 0.f) {
            alpha = alpha_set; //EM like when Posteriori estimation is null
          } else {
            alpha = 0.f; //Wiener like punctual supression when Posteriori estimation is not null
          }

          if(*(prev_frame) == 1) {
            Rprio = (1.f-alpha)*Rpost + alpha*gain_prev[k]*gain_prev[k]*(p2_prev[k]/noise_spectrum[k]);
          }else{
            Rprio = Rpost;
          }

          if(p2[k] > masked[k]) {
            gain = MAX(masked[k]/p2[k], Rprio/(Rprio+1.f)); //Rprio/(Rprio+1.0)
          } else {
            gain = 1.f; //Is signal worth of keeping
          }

          p2_prev[k] = p2[k];
          gain_prev[k] = gain;
          *(prev_frame) = 1;
          break;
      }

      //Apply over sustraction - This is intended to be used with powerspectrum
      //sustraction but it works with Wiener too

      //Limit gain to be applied
      Fk = over_reduc*(1.f-gain);

      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      Gk[k] =  1.f - Fk;
      Gk_prev[k] = Gk[k];

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for

}
