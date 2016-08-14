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

// //Circular Convolution
// void conv(float* x, float* h, float* y, int m, int n){
//   int i,j,k;
//   float a[m],x2[m];
//   if(m-n!=0){/*If length of both sequences are not equal*/
//     if(m>n) {/* Pad the smaller sequence with zero*/
//       for(i=n;i<m;i++)
//         h[i]=0;
//
//       n=m;
//     }
//     for(i=m;i<n;i++)
//       x[i]=0;
//
//     m=n;
//   }
//
//   y[0]=0;
//   a[0]=h[0];
//
//   for(j=1;j<n;j++) /*folding h(n) to h(-n)*/
//     a[j]=h[n-j];
//
//   /*Circular convolution*/
//   for(i=0;i<n;i++)
//     y[0]+=x[i]*a[i];
//
//   for(k=1;k<n;k++){
//     y[k]=0;
//     /*circular shift*/
//     for(j=1;j<n;j++)
//       x2[j]=a[j-1];
//
//     x2[0]=a[n-1];
//
//     for(i=0;i<n;i++){
//       a[i]=x2[i];
//       y[k]+=x[i]*x2[i];
//     }
//   }
// }
//
//
// void post_filter (float scale_factor,float* p2, float* Gk, float SNR_thresh, int fft_size_2) {
//   float SNR_estim, num = 0.f, den = 0.f, N, pf_thresh;
//   float post_filter[fft_size_2+1], Gk_copy[fft_size_2+1];
//   int k;
//
//   //Estimate Post-Filter
//   for (k = 0; k <= fft_size_2 ; k++) {
//     num += powf(fabs(Gk[k]*p2[k]),2);
//     den += powf(fabs(p2[k]),2);
//   }
//
//   SNR_estim = num/(den+FLT_MIN);
//
//   if (SNR_estim >= SNR_thresh){
//     pf_thresh = 1.f;
//   }else{
//     pf_thresh = SNR_estim;
//   }
//
//   if (pf_thresh == 1) {
//     N = 1.f;
//   } else {
//     N = 2.f*(1.f+(roundf(pf_thresh/SNR_thresh)*scale_factor)) + 1.f;
//   }
//
//   for (k = 0; k <= fft_size_2 ; k++) {
//     if(k<N){
//       post_filter[k] = 1.f/N;
//     }else{
//       post_filter[k] = 0.f;
//     }
//     Gk_copy[k] = fabs(Gk[k]);
//   }
//
//   //Apply Post-Filter
//   conv(Gk_copy,post_filter,Gk,fft_size_2+1,fft_size_2+1);
// }

void denoise_gain_ps(float over_reduc,
                  float* p2,
                  int fft_size_2,
                  float* Gk,
                  float* noise_thresholds) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (noise_thresholds[k] > FLT_MIN){
      // Power Subtraction
      if(p2[k] > FLT_MIN)
        gain = MAX(p2[k]-noise_thresholds[k], 0.f)/p2[k];
      else
        gain = 0.f;

      //Apply over sustraction
      Fk = over_reduc*(1.f-gain);

      //Limit gain to be applied
      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      gain =  1.f - Fk;

      //Assing gain to gain spectrum
      Gk[k] =  gain;

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for
}
void denoise_gain_w(float over_reduc,
                  float* p2,
                  int fft_size_2,
                  float* Gk,
                  float* noise_thresholds) {
  int k;
  float gain, Fk;

  for (k = 0; k <= fft_size_2 ; k++) {
    if (noise_thresholds[k] > FLT_MIN){
      //wiener estimation
      float rest = p2[k] - noise_thresholds[k];

      if(p2[k] > noise_thresholds[k])
        gain = (p2[k]) / (rest+noise_thresholds[k]);
      else
        gain = 0.f;

      //Apply over sustraction
      Fk = over_reduc*(1.f-gain);

      //Limit gain to be applied
      if(Fk < 0.f) Fk = 0.f;
      if(Fk > 1.f) Fk = 1.f;

      gain =  1.f - Fk;

      //Assing gain to gain spectrum
      Gk[k] =  gain;

    } else {
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  } //for
}
void denoise_gain_em(float over_reduc,
                        float alpha_set,
                        int* prev_frame,
                        float* p2,
                        float* p2_prev,
                        float* gain_prev,
                        int fft_size_2,
                        float* Gk,
                        float* noise_thresholds) {
    int k;
    float gain, Fk, Rpost, Rprio;

    for (k = 0; k <= fft_size_2 ; k++) {
      if (noise_thresholds[k] > FLT_MIN){
        // EM using Wolfe and Godsill optimization
        Rpost = MAX(p2[k]/noise_thresholds[k]-1.f, 0.f);

        if(*(prev_frame) == 1) {
          Rprio = (1.f-alpha_set)*Rpost + alpha_set*gain_prev[k]*gain_prev[k]*(p2_prev[k]/noise_thresholds[k]);
        }else{
          Rprio = Rpost;
        }

        //Ephraim-Malah noise suppression, from Godsill and Wolfe 2001 paper (cheaper)
        float r = MAX(Rprio/(1.f+Rprio),FLT_MIN);
        float V = (Rprio/(1.f+Rprio))*(Rpost+1.f);
        gain = sqrtf( r * (1.f+V)/(Rpost+1.f) );

        p2_prev[k] = p2[k];
        gain_prev[k] = gain;
        *(prev_frame) = 1;

        //Apply over sustraction
        Fk = over_reduc*(1.f-gain);

        //Limit gain to be applied
        if(Fk < 0.f) Fk = 0.f;
        if(Fk > 1.f) Fk = 1.f;

        gain =  1.f - Fk;

        //Assing gain to gain spectrum
        Gk[k] =  gain;

      } else {
        //Otherwise we keep everything as is
        Gk[k] = 1.f;
      }
    } //for
  }
void denoise_gain_cmsr(float over_reduc,
                        float alpha_set,
                        int* prev_frame,
                        float* p2,
                        float* p2_prev,
                        float* gain_prev,
                        int fft_size_2,
                        float* Gk,
                        float* noise_thresholds) {
    int k;
    float gain, Fk, Rpost, Rprio, alpha;

    for (k = 0; k <= fft_size_2 ; k++) {
      if (noise_thresholds[k] > FLT_MIN){

        // CMSR (modified EM)
        Rpost = MAX(p2[k]/noise_thresholds[k]-1.f, 0.f);

        if (Rpost > 0.f) {
          alpha = alpha_set; //EM like when Posteriori estimation is null
        } else {
          alpha = 0.f; //Wiener like punctual supression when Posteriori estimation is not null
        }

        if(*(prev_frame) == 1) {
          Rprio = (1.f-alpha)*Rpost + alpha*gain_prev[k]*gain_prev[k]*(p2_prev[k]/noise_thresholds[k]);
        }else{
          Rprio = Rpost;
        }

        //Ephraim-Malah noise suppression, from Godsill and Wolfe 2001 paper (cheaper)
        float r = MAX(Rprio/(1.f+Rprio),FLT_MIN);
        float V = (Rprio/(1.f+Rprio))*(Rpost+1.f);
        gain = sqrtf( r * (1.f+V)/(Rpost+1.f) );

        p2_prev[k] = p2[k];
        gain_prev[k] = gain;
        *(prev_frame) = 1;

        //Apply over sustraction
        Fk = over_reduc*(1.f-gain);

        //Limit gain to be applied
        if(Fk < 0.f) Fk = 0.f;
        if(Fk > 1.f) Fk = 1.f;

        gain =  1.f - Fk;

        //Assing gain to gain spectrum
        Gk[k] =  gain;

      } else {
        //Otherwise we keep everything as is
        Gk[k] = 1.f;
      }
    } //for
  }
