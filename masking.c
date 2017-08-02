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

//extra values
#define N_BARK_BANDS 25
#define USE_POWER_SPECTRUM_FOR_TONALITY false
#define EXPONENT 10
#define HIGH_FREQ_BIAS 7.f

#define MASKING 6.0
#define GAIN_SMOOTH 0.5f //smoothing of gain proposed in virag method

#define ARRAYACCESS(a, i, j) ((a)[(i) * N_BARK_BANDS + (j)]) //This is for SSF Matrix recall

//Proposed by Sinha and Tewfik and explained by Virag
const float relative_thresholds[N_BARK_BANDS] = {
  -16.f,
  -17.f,
  -18.f,
  -19.f,
  -20.f,
  -21.f,
  -22.f,
  -23.f,
  -24.f,
  -25.f,
  -25.f,
  -25.f,
  -25.f,
  -25.f,
  -25.f,
  -24.f,
  -23.f,
  -22.f,
  -19.f,
  -18.f,
  -18.f,
  -18.f,
  -18.f,
  -18.f,
  -18.f,
};

//fft to bark bilinear transform
void compute_bark_mapping(float* bark_z,int fft_size_2, int srate) {
  int k;
  float freq;
  /* compute the bark z value for this frequency bin */
  for(k = 0 ; k <= fft_size_2 ; k++) {
    freq = (float)srate / 2.f /(float)(fft_size_2)*(float)k ; //bin to freq
    //bark_z[k] = 7.f*logf(freq/650.f + sqrtf(1.f + (freq/650.f)*(freq/650.f))) ;
    bark_z[k] = 1.f + 13.f*atanf(0.00076f*freq) + 3.5f*atanf(powf(freq/7500,2)) ;
  }
}

void compute_absolute_thresholds(float* absolute_thresholds,int fft_size_2, int srate) {
  int k;
  float freq;
  for(k = 1 ; k <= fft_size_2 ; k++) {
    freq = Index2Freq(k,srate,fft_size_2); //bin to freq
    absolute_thresholds[k] = 3.64*powf((freq/1000),-0.8) - 6.5*exp(-0.6*powf((freq/1000 - 3.3),2)) + powf(10,-3)*powf((freq /1000),4);//dB scale
  }
}

void compute_SSF(float* SSF) {
  int i,j;
  float y;
  for(i = 0 ; i < N_BARK_BANDS ; i++) {
    for(j = 0 ; j < N_BARK_BANDS ; j++) {
      y = (i+1)-(j+1);
      //Spreading function (Schroeder)
      ARRAYACCESS(SSF,i,j) = 15.81 + 7.5*(y +0.474) - 17.5*sqrtf(1.f+(y +0.474)*(y +0.474));//dB scale
      //db to Linear
      ARRAYACCESS(SSF,i,j) = powf(10.f,ARRAYACCESS(SSF,i,j)/10.f);
    }
  }
}

//Convolution by multiplication of a Toepliz matrix to a vector
void convolve_with_SSF(float* SSF,//signal
float* bark_spectrum,//kernel
float* spreaded_spectrum){

  int i,j;
  for (i = 0; i < N_BARK_BANDS; i++){
    spreaded_spectrum[i] = 0.f;
    for (j = 0; j < N_BARK_BANDS; j++){
      spreaded_spectrum[i] += ARRAYACCESS(SSF,i,j)*bark_spectrum[j];
    }
  }
}

//Computes the energy of each bark band
void compute_bark_spectrum(float* bark_z,
float* bark_spectrum,
float* spectrum,
int* intermediate_band_bins,
int* n_bins_per_band){
  int j;
  int last_position = 0;

  //Critical band Analysis
  for (j = 0; j < N_BARK_BANDS; j++){
   //Using mapping to bark previously computed

   int cont = 0;
   if(j==0) cont = 1; //Do not take into account the DC component

   bark_spectrum[j] = 0.f;
   //If we are on the same band for the bin
   while(floor(bark_z[last_position+cont]) == (j+1)){//First bark band is 1
     bark_spectrum[j] += spectrum[last_position+cont];
     cont++;
   }

   //Dividing the energy in the bark band with the number of bins of the band
   // bark_spectrum[j] /= cont;

   //Move the position to the next group of bins from the upper bark band
   last_position += cont;

   //store bin information
   n_bins_per_band[j] = cont;
   intermediate_band_bins[j] = last_position;
  }
}

//Computes the tonality factor using the spectral flatness
float compute_tonality_factor(float* bark_spectrum,
float* spectrum,
float fft_size_2,
int* n_bins_per_band){

  int j, k;
  float SFM, SFM_array[N_BARK_BANDS],Gm[N_BARK_BANDS],Am[N_BARK_BANDS],tonality_factor;
  float sum_p = 0.f,sum_log_p = 0.f;

  if(USE_POWER_SPECTRUM_FOR_TONALITY){
    //Using power spectrum to compute the tonality factor
    for (k = 0; k <= fft_size_2; k++) {
      //For spectral flatness measures
      sum_p += spectrum[k];
      sum_log_p += log10f(spectrum[k]);
    }
    //spectral flatness measure using Geometric and Arithmetic means of the spectrum cleaned previously
    //Using log propieties and definition of spectral flatness https://en.wikipedia.org/wiki/Spectral_flatness
    //Robinson thesis explains this reexpresion in detail
    SFM = 10.f*(sum_log_p/(float)fft_size_2 - log10f(sum_p/(float)fft_size_2));//this value is in db scale

    //Tonality factor in db scale
    tonality_factor = MIN(SFM/-60.f, 1.f);

  }else{
    //Using bark spectrum to compute tonality factor
    for (j = 0; j < N_BARK_BANDS; j++){
      //For spectral flatness measures (Geometric and Arithmetic mean)
      Gm[j] = powf(bark_spectrum[j],1.f/(float)(n_bins_per_band[j]+1));
      Am[j] = 1.f/(float)(n_bins_per_band[j]+1)*bark_spectrum[j];

      //spectral flatness measure using Geometric and Arithmetic means
      SFM_array[j] = 10.f*log10f(Gm[j]/Am[j]);//this value is in db scale

      //taking into account max db SPL
      SFM_array[j] /= -60.f;
    }

    //Tonality factor in db scale
    tonality_factor = MIN(min_spectral_value(SFM_array,N_BARK_BANDS), 1.f);
  }

  return tonality_factor;
}

//masking threshold calculation
void compute_masking_thresholds(float* bark_z,
float* absolute_thresholds,
float* SSF,
float* spectrum,
int fft_size_2,
float* masking_thresholds) {
  int k, j, start, end;
  int intermediate_band_bins[N_BARK_BANDS];
  int n_bins_per_band[N_BARK_BANDS];
  float bark_spectrum[N_BARK_BANDS];
  float threshold_j[N_BARK_BANDS];
  float masking_offset[N_BARK_BANDS];
  float spreaded_spectrum[N_BARK_BANDS];
  float unity_gain_bark_spectrum[N_BARK_BANDS];
  memset(unity_gain_bark_spectrum, 1.f,N_BARK_BANDS); //initialization
  float spreaded_unity_gain_bark_spectrum[N_BARK_BANDS];
  float tonality_factor;

  //First we get the energy in each bark band
  compute_bark_spectrum(bark_z,
                        bark_spectrum,
                        spectrum,
                        intermediate_band_bins,
                        n_bins_per_band);

  //Now that we have the bark spectrum

  //Then we compute the tonality_factor for each band
  tonality_factor = compute_tonality_factor(bark_spectrum,
                                            spectrum,
                                            fft_size_2,
                                            n_bins_per_band);

  //Convolution bewtween the bark spectrum and SSF (Toepliz matrix multiplication)
  convolve_with_SSF(SSF,bark_spectrum,spreaded_spectrum);
  //Convolve unitary energy bark spectrum with SSF
  convolve_with_SSF(SSF,unity_gain_bark_spectrum,spreaded_unity_gain_bark_spectrum);

  for (j = 0; j < N_BARK_BANDS; j++){

    //Masking offset
    masking_offset[j] = (tonality_factor*(14.5+(j+1))+5.5*(1.f - tonality_factor));
    //masking_offset[j] = relative_thresholds[j];

    //Consider tonal noise in upper bands (j>15) due to musical noise of the power Sustraction
    //that was used at First
    //if(j>15) masking_offset[j] -= HIGH_FREQ_BIAS;

    //spread Masking threshold
    threshold_j[j] = powf(10.f,log10f(spreaded_spectrum[j]) -  masking_offset[j]/10.f);
    //threshold_j[j] = 10.f*log10f(spreaded_spectrum[j]) +  masking_offset[j];

    //Renormalization
    threshold_j[j] -= 10.f*log10f(spreaded_unity_gain_bark_spectrum[j]);
    //threshold_j[j] *= (1.f/(10*log10f(spreaded_unity_gain_bark_spectrum[j])));
    //printf("%f\n",threshold_j[j]);

    //Relating the spread masking threshold to the critical band masking thresholds

    //Border case
    if(j == 0){
      start = 0;
    } else {
      start = intermediate_band_bins[j-1];
    }
    end = intermediate_band_bins[j];

    //printf("%f\n",spreaded_spectrum[j]);

    for(k = start; k < end; k++){
      //Taking into account the absolute hearing thresholds (90 dbfs of reference level)
      masking_thresholds[k] = MAX(threshold_j[j],absolute_thresholds[k]);
      //masking_thresholds[k] = threshold_j[j] + (90.f - absolute_thresholds[k]);
      //masking_thresholds[k] = threshold_j[j];
    }
  }
}

//alpha and beta computation to be used in general spectral Sustraction
void compute_alpha_and_beta(float* fft_p2,
float* noise_thresholds_p2,
int fft_size_2,
float* alpha,
//float* beta,
float* bark_z,
float* absolute_thresholds,
float* SSF,
float* max_masked,
float* min_masked,
float masking) {

  int k;
  float masking_thresholds[fft_size_2+1];
  float estimated_clean[fft_size_2+1];


  //Noise masking threshold must be computed from a clean signal
  //therefor we aproximate a clean signal using a power Sustraction over
  //the original noisy one

  //basic Power Sustraction to estimate clean signal
  for (k = 0; k <= fft_size_2; k++) {
    estimated_clean[k] =  MAX(fft_p2[k]-noise_thresholds_p2[k],FLT_MIN);
  }

  //Now we can compute noise masking threshold from this clean signal
  compute_masking_thresholds(bark_z,
                             absolute_thresholds,
                             SSF,
                             estimated_clean,
                             fft_size_2,
                             masking_thresholds);

  /*Get alpha and beta based on masking thresholds
  *beta and alpha values would adapt based on masking thresholds
  *frame to frame for optimal oversustraction and noise floor parameter in each one
  *noise floor would better be controled by user using the amount of reduction
  *so beta is not modified
  */

  //First we need the maximun and the minimun value of the masking threshold
  *(max_masked) = MAX(max_spectral_value(masking_thresholds,fft_size_2),*(max_masked));
  *(min_masked) = MIN(min_spectral_value(masking_thresholds,fft_size_2),*(min_masked));

  printf("%f\n",*(max_masked) );
  printf("%f\n",*(min_masked) );

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
       //alpha[k] = powf(1 - (masking_thresholds[k] - ALPHA_MIN)/(masking - ALPHA_MIN),EXPONENT) * (masking - ALPHA_MIN) + masking;
       alpha[k] = ALPHA_MIN + (masking - ALPHA_MIN)/(*(min_masked) - *(max_masked)) * (masking_thresholds[k]- ALPHA_MIN);
       //beta[k] = (1 - (masking_thresholds[k] - BETA_MIN)/(masking - BETA_MIN)) * (masking - BETA_MIN) + masking;
    }
  }
}
