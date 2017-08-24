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
#define HIGH_FREQ_BIAS 20.f
#define S_AMP 1.f
#define AT_SINE_WAVE_FREQ 1000.f
#define REFERENCE_LEVEL 90.f //dbSPL level of reproduction

#define ARRAYACCESS(a, i, j) ((a)[(i) * N_BARK_BANDS + (j)]) //This is for SSF Matrix recall

//Proposed by Sinha and Tewfik and explained by Virag
const float relative_thresholds[N_BARK_BANDS] = { -16.f, -17.f, -18.f, -19.f, -20.f, -21.f, -22.f, -23.f, -24.f, -25.f, -25.f, -25.f, -25.f, -25.f, -25.f, -24.f, -23.f, -22.f, -19.f, -18.f, -18.f, -18.f, -18.f, -18.f, -18.f};

//fft to bark bilinear transform
void
compute_bark_mapping(float* bark_z,int fft_size_2, int srate)
{
  int k;
  float freq;
  /* compute the bark z value for this frequency bin */
  for(k = 0 ; k <= fft_size_2 ; k++)
  {
    freq = (float)srate / 2.f /(float)(fft_size_2)*(float)k ; //bin to freq
    //bark_z[k] = 7.f*logf(freq/650.f + sqrtf(1.f + (freq/650.f)*(freq/650.f))) ;
    bark_z[k] = 1.f + 13.f*atanf(0.00076f*freq) + 3.5f*atanf(powf(freq/7500.f,2.f)) ;
  }
}

void
compute_SSF(float* SSF)
{
  int i,j;
  float y;
  for(i = 0 ; i < N_BARK_BANDS ; i++)
  {
    for(j = 0 ; j < N_BARK_BANDS ; j++)
    {
      y = (i+1)-(j+1);
      //Spreading function (Schroeder)
      ARRAYACCESS(SSF,i,j) = 15.81 + 7.5*(y +0.474) - 17.5*sqrtf(1.f+(y +0.474)*(y +0.474));//dB scale
      //db to Linear
      ARRAYACCESS(SSF,i,j) = powf(10.f,ARRAYACCESS(SSF,i,j)/10.f);
    }
  }
}

//Convolution by multiplication of a Toepliz matrix to a vector
void convolve_with_SSF(float* SSF, float* bark_spectrum, float* spreaded_spectrum)
{
  int i,j;
  for (i = 0; i < N_BARK_BANDS; i++)
  {
    spreaded_spectrum[i] = 0.f;
    for (j = 0; j < N_BARK_BANDS; j++)
    {
      spreaded_spectrum[i] += ARRAYACCESS(SSF,i,j)*bark_spectrum[j];
    }
  }
}

//Computes the energy of each bark band
void
compute_bark_spectrum(float* bark_z, float* bark_spectrum, float* spectrum,
                      float* intermediate_band_bins, float* n_bins_per_band)
{
  int j;
  int last_position = 0;

  //Critical band Analysis
  for (j = 0; j < N_BARK_BANDS; j++)
  {
   //Using mapping to bark previously computed

   int cont = 0;
   if(j==0) cont = 1; //Do not take into account the DC component

   bark_spectrum[j] = 0.f;
   //If we are on the same band for the bin
   while(floor(bark_z[last_position+cont]) == (j+1))
   {//First bark band is 1
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

void
spl_reference(float* spl_reference_values, int fft_size_2, int srate,
              float* input_fft_buffer_at, float* output_fft_buffer_at,
              fftwf_plan* forward_at)
{
  int k;
  float sinewave[2*fft_size_2];
  float window[2*fft_size_2];
  float fft_p2_at[fft_size_2+1];
  float fft_magnitude_at[fft_size_2+1];
  float fft_phase_at[fft_size_2+1];
  float fft_p2_at_dbspl[fft_size_2+1];

  //Generate a fullscale sine wave of 1 kHz
  for(k = 0 ; k < 2*fft_size_2 ; k++)
  {
    sinewave[k] = S_AMP * sinf((2.f * M_PI * k * AT_SINE_WAVE_FREQ) / (float)srate);
  }

  //Windowing the sinewave
  fft_window(window, 2*fft_size_2,0);//von-Hann window
  for(k = 0 ; k < 2*fft_size_2 ; k++)
  {
    input_fft_buffer_at[k] = sinewave[k] * window[k];
  }

  //Do FFT
  fftwf_execute(*forward_at);

  //Get the magnitude
  get_info_from_bins(fft_p2_at, fft_magnitude_at, fft_phase_at, fft_size_2, 2*fft_size_2,
                     output_fft_buffer_at);

  //Convert to db and taking into account 90dbfs of reproduction loudness
  for(k = 0 ; k <= fft_size_2 ; k++)
  {
    fft_p2_at_dbspl[k] = REFERENCE_LEVEL - 10.f*log10f(fft_p2_at[k]);
  }

  memcpy(spl_reference_values,fft_p2_at_dbspl,sizeof(float)*(fft_size_2+1));

}

void
convert_to_dbspl(float* spl_reference_values,float* masking_thresholds, int fft_size_2)
{
  for(int k = 0 ; k <= fft_size_2 ; k++)
  {
    masking_thresholds[k] += spl_reference_values[k];
  }
}

void
compute_absolute_thresholds(float* absolute_thresholds,int fft_size_2, int srate)
{
  int k;
  float freq;

  for(k = 1 ; k <= fft_size_2 ; k++)
  {//As explained by thiemann
    freq = bin_to_freq(k,srate,fft_size_2); //bin to freq
    absolute_thresholds[k] = 3.64*powf((freq/1000),-0.8) - 6.5*exp(-0.6*powf((freq/1000 - 3.3),2)) + powf(10,-3)*powf((freq /1000),4);//dBSPL scale
  }
}

//Computes the tonality factor using the spectral flatness
float
compute_tonality_factor(float* spectrum, float* intermediate_band_bins,
                        float* n_bins_per_band, int band)
{
  int k;
  float SFM, tonality_factor;
  float sum_p = 0.f,sum_log_p = 0.f;
  int start_pos, end_pos = 0;

  //Mapping to bark bands
  if(band == 0)
  {
    start_pos = band;
    end_pos = n_bins_per_band[band];
  }
  else
  {
    start_pos = intermediate_band_bins[band-1];
    end_pos = intermediate_band_bins[band-1] + n_bins_per_band[band];
  }

  //Using power spectrum to compute the tonality factor
  for (k = start_pos; k < end_pos; k++)
  {
    //For spectral flatness measures
    sum_p += spectrum[k];
    sum_log_p += log10f(spectrum[k]);
  }
  //spectral flatness measure using Geometric and Arithmetic means of the spectrum cleaned previously
  //Using log propieties and definition of spectral flatness https://en.wikipedia.org/wiki/Spectral_flatness
  //Robinson thesis explains this reexpresion in detail
  SFM = 10.f*(sum_log_p/(float)(n_bins_per_band[band]) - log10f(sum_p/(float)(n_bins_per_band[band])));//this value is in db scale

  //Tonality factor in db scale
  tonality_factor = MIN(SFM/-60.f, 1.f);

  return tonality_factor;
}

//masking threshold calculation
void
compute_masking_thresholds(float* bark_z, float* absolute_thresholds, float* SSF,
                           float* spectrum, int fft_size_2, float* masking_thresholds,
                           float* spreaded_unity_gain_bark_spectrum,
                           float* spl_reference_values)
{
  int k, j, start_pos, end_pos;
  float intermediate_band_bins[N_BARK_BANDS];
  float n_bins_per_band[N_BARK_BANDS];
  float bark_spectrum[N_BARK_BANDS];
  float threshold_j[N_BARK_BANDS];
  float masking_offset[N_BARK_BANDS];
  float spreaded_spectrum[N_BARK_BANDS];
  float tonality_factor;

  //First we get the energy in each bark band
  compute_bark_spectrum(bark_z, bark_spectrum, spectrum, intermediate_band_bins,
                        n_bins_per_band);

  //Now that we have the bark spectrum
  //Convolution bewtween the bark spectrum and SSF (Toepliz matrix multiplication)
  convolve_with_SSF(SSF, bark_spectrum, spreaded_spectrum);

  for (j = 0; j < N_BARK_BANDS; j++)
  {
    //Then we compute the tonality_factor for each band (1 tone like 0 noise like)
    tonality_factor = compute_tonality_factor(spectrum, intermediate_band_bins,  n_bins_per_band, j);//Uses power spectrum

    //Masking offset
    masking_offset[j] = (tonality_factor*(14.5+(j+1)) + 5.5*(1.f - tonality_factor));

    //Using offset proposed by Virag
    //masking_offset[j] = relative_thresholds[j];

    //Consider tonal noise in upper bands (j>15) due to musical noise of the power Sustraction that was used at First
    //if(j>15) masking_offset[j] += HIGH_FREQ_BIAS;

    //spread Masking threshold
    threshold_j[j] = powf(10.f,log10f(spreaded_spectrum[j]) - (masking_offset[j]/10.f));

    //Renormalization
    threshold_j[j] -= 10.f*log10f(spreaded_unity_gain_bark_spectrum[j]);

    //Relating the spread masking threshold to the critical band masking thresholds
    //Border case
    if(j == 0)
    {
      start_pos = 0;
    }
    else
    {
      start_pos = intermediate_band_bins[j-1];
    }
    end_pos = intermediate_band_bins[j];

    for(k = start_pos; k < end_pos; k++)
    {
      masking_thresholds[k] = threshold_j[j];
    }
  }

  //Masking thresholds need to be converted to db spl scale in order to be compared with
  //absolute threshold of hearing
  convert_to_dbspl(spl_reference_values,masking_thresholds,fft_size_2);

  //Take into account the absolute_thresholds of hearing
  for(k = 0; k <= fft_size_2; k++)
  {
    masking_thresholds[k] = MAX(masking_thresholds[k],absolute_thresholds[k]);
  }
}

//alpha and beta computation according to Virag
void
compute_alpha_and_beta(float* fft_p2, float* noise_thresholds_p2, int fft_size_2,
                       float* alpha_masking, float* beta_masking, float* bark_z,
                       float* absolute_thresholds, float* SSF,
                       float* spreaded_unity_gain_bark_spectrum,
                       float* spl_reference_values, float masking_value,
                       float reduction_value)
{
  int k;
  float masking_thresholds[fft_size_2+1];
  float estimated_clean[fft_size_2+1];
  float normalized_value;


  //Noise masking threshold must be computed from a clean signal
  //therefor we aproximate a clean signal using a power Sustraction over
  //the original noisy one

  //basic Power Sustraction to estimate clean signal
  for (k = 0; k <= fft_size_2; k++)
  {
    estimated_clean[k] =  MAX(fft_p2[k] - noise_thresholds_p2[k],FLT_MIN);
  }

  //Now we can compute noise masking threshold from this clean signal
  compute_masking_thresholds(bark_z, absolute_thresholds, SSF, estimated_clean,
                             fft_size_2, masking_thresholds,
                             spreaded_unity_gain_bark_spectrum, spl_reference_values);

  /*Get alpha and beta based on masking thresholds
  *beta and alpha values would adapt based on masking thresholds
  *frame to frame for optimal oversustraction and noise floor parameter in each one
  *noise floor would better be controled by user using the amount of reduction
  *so beta is not modified
  */

  //First we need the maximun and the minimun value of the masking threshold
  float max_masked_tmp = max_spectral_value(masking_thresholds,fft_size_2);
  float min_masked_tmp = min_spectral_value(masking_thresholds,fft_size_2);

  for (k = 0; k <= fft_size_2; k++)
  {
    //new alpha and beta vector
    if(masking_thresholds[k] == max_masked_tmp)
    {
       alpha_masking[k] = ALPHA_MIN;
       beta_masking[k] = BETA_MIN;
    }
    if(masking_thresholds[k] == min_masked_tmp)
    {
       alpha_masking[k] = masking_value;
       beta_masking[k] = reduction_value;
    }
    if(masking_thresholds[k] < max_masked_tmp && masking_thresholds[k] > min_masked_tmp)
    {
      //Linear interpolation of the value between max and min masked threshold values
      normalized_value = (masking_thresholds[k]-min_masked_tmp)/(max_masked_tmp-min_masked_tmp);

      alpha_masking[k] = (1.f - normalized_value)*ALPHA_MIN + normalized_value*masking_value;
      beta_masking[k] = (1.f - normalized_value)*BETA_MIN + normalized_value*reduction_value;
    }
  }
}
