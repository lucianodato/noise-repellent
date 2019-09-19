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

/**
* \file masking.c
* \author Luciano Dato
* \brief Methods for masking threshold estimation
*/

#include <float.h>
#include <math.h>
#include <fftw3.h>

//masking thresholds values recomended by virag
#define ALPHA_MAX 6.f
#define ALPHA_MIN 1.f
#define BETA_MAX 0.02f
#define BETA_MIN 0.0f

//extra values
#define N_BARK_BANDS 25
#define AT_SINE_WAVE_FREQ 1000.f
#define REFERENCE_LEVEL 90.f //dbSPL level of reproduction

#define BIAS 0
#define HIGH_FREQ_BIAS 20.f
#define S_AMP 1.f

#define ARRAYACCESS(a, i, j) ((a)[(i)*N_BARK_BANDS + (j)]) //This is for SSF Matrix recall

//Proposed by Sinha and Tewfik and explained by Virag
static const float
relative_thresholds[N_BARK_BANDS] = {-16.f, -17.f, -18.f, -19.f, -20.f, -21.f, -22.f, -23.f, -24.f, -25.f, -25.f, -25.f, -25.f, -25.f, -25.f, -24.f, -23.f, -22.f, -19.f, -18.f, -18.f, -18.f, -18.f, -18.f, -18.f};

/**
* Fft to bark bilinear scale transform. This computes the corresponding bark band for
* each fft bin and generates an array that
* inicates this mapping.
* \param bark_z defines the bark to linear mapping for current spectrum config
* \param fft_size_2 is half of the fft size
* \param srate current sample rate of the host
*/
static void
compute_bark_mapping(float *bark_z, int fft_size_2, int srate)
{
  int k;
  float freq;

  for (k = 0; k <= fft_size_2; k++)
  {
    freq = (float)srate / (2.f * (float)(fft_size_2) * (float)k); //bin to freq
    bark_z[k] = 1.f + 13.f * atanf(0.00076f * freq) + 3.5f * atanf(powf(freq / 7500.f, 2.f));
  }
}

/**
* Computes the spectral spreading function of Schroeder as a matrix using bark scale.
* This is to perform a convolution between this function and a bark spectrum. The
* complete explanation for this is in Robinsons master thesis 'Perceptual model for
* assessment of coded audio'.
* \param SSF defines the spreading function matrix
*/
static void
compute_SSF(float *SSF)
{
  int i, j;
  float y;
  for (i = 0; i < N_BARK_BANDS; i++)
  {
    for (j = 0; j < N_BARK_BANDS; j++)
    {
      y = (i + 1) - (j + 1);
      //Spreading function (Schroeder)
      ARRAYACCESS(SSF, i, j) = 15.81f + 7.5f * (y + 0.474f) - 17.5f * sqrtf(1.f + (y + 0.474f) * (y + 0.474f)); //dB scale
      //db to Linear
      ARRAYACCESS(SSF, i, j) = powf(10.f, ARRAYACCESS(SSF, i, j) / 10.f);
    }
  }
}

/**
* Convolution between the spreading function by multiplication of a Toepliz matrix
* to a bark spectrum. The complete explanation for this is in Robinsons master thesis
* 'Perceptual model for assessment of coded audio'.
* \param SSF defines the spreading function matrix
* \param bark_spectrum the bark spectrum values of current power spectrum
* \param spreaded_spectrum result of the convolution bewtween SSF and the bark spectrum
*/
static void
convolve_with_SSF(float *SSF, float *bark_spectrum, float *spreaded_spectrum)
{
  int i, j;
  for (i = 0; i < N_BARK_BANDS; i++)
  {
    spreaded_spectrum[i] = 0.f;
    for (j = 0; j < N_BARK_BANDS; j++)
    {
      spreaded_spectrum[i] += ARRAYACCESS(SSF, i, j) * bark_spectrum[j];
    }
  }
}

/**
* Computes the energy of each bark band taking a power or magnitude spectrum. It performs
* the mapping from linear fft scale to the bark scale. Often called critical band Analysis
* \param bark_z defines the bark to linear mapping for current spectrum config
* \param bark_spectrum the bark spectrum values of current power spectrum
* \param spectrum is the power spectum array
* \param intermediate_band_bins holds the bin numbers that are limits of each band
* \param n_bins_per_band holds the the number of bins in each band
*/
static void
compute_bark_spectrum(float *bark_z, float *bark_spectrum, float *spectrum,
                      float *intermediate_band_bins, float *n_bins_per_band)
{
  int j;
  int last_position = 0;

  for (j = 0; j < N_BARK_BANDS; j++)
  {
    int cont = 0;
    if (j == 0)
      cont = 1; //Do not take into account the DC component

    bark_spectrum[j] = 0.f;
    //If we are on the same band for the bin
    while (floor(bark_z[last_position + cont]) == (j + 1))
    { //First bark band is 1
      bark_spectrum[j] += spectrum[last_position + cont];
      cont++;
    }
    //Move the position to the next group of bins from the upper bark band
    last_position += cont;

    //store bin information
    n_bins_per_band[j] = cont;
    intermediate_band_bins[j] = last_position;
  }
}

/**
* Computes the reference spectrum to perform the db to dbSPL conversion. It uses
* a full scale 1khz sine wave as the reference signal and performs an fft transform over
* it to get the magnitude spectrum and then scales it using a reference dbSPL level. This
* is used because the absolute thresholds of hearing are obtained in SPL scale so it's
* necessary to compare this thresholds with obtained masking thresholds using the same
* scale.
* \param spl_reference_values defines the reference values for each bin to convert from db to db SPL
* \param fft_size_2 is half of the fft size
* \param srate current sample rate of the host
* \param input_fft_buffer_at input buffer for the reference sinewave fft transform
* \param output_fft_buffer_at output buffer for the reference sinewave fft transform
* \param forward_at fftw plan for the reference sinewave fft transform
*/
static void
spl_reference(float *spl_reference_values, int fft_size_2, int srate,
              float *input_fft_buffer_at, float *output_fft_buffer_at,
              fftwf_plan *forward_at)
{
  int k;
  float sinewave[2 * fft_size_2];
  float window[2 * fft_size_2];
  float fft_p2_at[fft_size_2 + 1];
  float fft_magnitude_at[fft_size_2 + 1];
  float fft_phase_at[fft_size_2 + 1];
  float fft_p2_at_dbspl[fft_size_2 + 1];

  //Generate a fullscale sine wave of 1 kHz
  for (k = 0; k < 2 * fft_size_2; k++)
  {
    sinewave[k] = S_AMP * sinf((2.f * M_PI * k * AT_SINE_WAVE_FREQ) / (float)srate);
  }

  //Windowing the sinewave
  fft_window(window, 2 * fft_size_2, 0); //von-Hann window
  for (k = 0; k < 2 * fft_size_2; k++)
  {
    input_fft_buffer_at[k] = sinewave[k] * window[k];
  }

  //Do FFT
  fftwf_execute(*forward_at);

  //Get the magnitude
  get_info_from_bins(fft_p2_at, fft_magnitude_at, fft_phase_at, fft_size_2, 2 * fft_size_2,
                     output_fft_buffer_at);

  //Convert to db and taking into account 90dbfs of reproduction loudness
  for (k = 0; k <= fft_size_2; k++)
  {
    fft_p2_at_dbspl[k] = REFERENCE_LEVEL - 10.f * log10f(fft_p2_at[k]);
  }

  memcpy(spl_reference_values, fft_p2_at_dbspl, sizeof(float) * (fft_size_2 + 1));
}

/**
* dB scale to dBSPL conversion. This is to convert masking thresholds from db to dbSPL
* scale to then compare them to absolute threshold of hearing.
* \param spl_reference_values defines the reference values for each bin to convert from db to db SPL
* \param masking_thresholds the masking thresholds obtained in db scale
* \param fft_size_2 is half of the fft size
*/
static void
convert_to_dbspl(float *spl_reference_values, float *masking_thresholds, int fft_size_2)
{
  for (int k = 0; k <= fft_size_2; k++)
  {
    masking_thresholds[k] += spl_reference_values[k];
  }
}

/**
* Computes the absolute thresholds of hearing to contrast with the masking thresholds.
* This formula is explained in Thiemann thesis 'Acoustic Noise Suppression for Speech
* Signals using Auditory Masking Effects'
* \param absolute_thresholds defines the absolute thresholds of hearing for current spectrum config
* \param fft_size_2 is half of the fft size
* \param srate current sample rate of the host
*/
static void
compute_absolute_thresholds(float *absolute_thresholds, int fft_size_2, int srate)
{
  int k;
  float freq;

  for (k = 1; k <= fft_size_2; k++)
  {                                                                                                                                                     //As explained by thiemann
    freq = bin_to_freq(k, srate, fft_size_2);                                                                                                           //bin to freq
    absolute_thresholds[k] = 3.64f * powf((freq / 1000.f), -0.8f) - 6.5f * exp(-0.6f * powf((freq / 1000.f - 3.3f), 2.f)) + powf(10.f, -3.f) * powf((freq / 1000.f), 4.f); //dBSPL scale
  }
}

/**
* Computes the tonality factor using the spectral flatness for a given bark band
* spectrum. Values are in db scale. An inferior limit of -60 db is imposed. To avoid zero
* logs some trickery is used as explained in https://en.wikipedia.org/wiki/Spectral_flatness
* Robinsons thesis explains this further too.
* \param spectrum is the power spectum array
* \param intermediate_band_bins holds the bin numbers that are limits of each band
* \param n_bins_per_band holds the the number of bins in each band
* \param band the bark band given
*/
static float
compute_tonality_factor(float *spectrum, float *intermediate_band_bins,
                        float *n_bins_per_band, int band)
{
  int k;
  float SFM, tonality_factor;
  float sum_p = 0.f, sum_log_p = 0.f;
  int start_pos, end_pos = 0;

  //Mapping to bark bands
  if (band == 0)
  {
    start_pos = band;
    end_pos = n_bins_per_band[band];
  }
  else
  {
    start_pos = intermediate_band_bins[band - 1];
    end_pos = intermediate_band_bins[band - 1] + n_bins_per_band[band];
  }

  //Using power spectrum to compute the tonality factor
  for (k = start_pos; k < end_pos; k++)
  {
    //For spectral flatness measures
    sum_p += spectrum[k];
    sum_log_p += log10f(spectrum[k]);
  }
  //spectral flatness measure using Geometric and Arithmetic means of the spectrum
  SFM = 10.f * (sum_log_p / (float)(n_bins_per_band[band]) - log10f(sum_p / (float)(n_bins_per_band[band]))); //this value is in db scale

  //Tonality factor in db scale
  tonality_factor = MIN(SFM / -60.f, 1.f);

  return tonality_factor;
}

/**
* Johnston Masking threshold calculation. This are greatly explained in Robinsons thesis
* and Thiemann thesis too and in Virags paper 'Single Channel Speech Enhancement Based on
* Masking Properties of the Human Auditory System'. Some optimizations suggested in
* Virags work are implemented but not used as they seem to not be necessary in modern
* age computers.
* \param bark_z defines the bark to linear mapping for current spectrum config
* \param absolute_thresholds defines the absolute thresholds of hearing for current spectrum config
* \param SSF defines the spreading function matrix
* \param spectrum is the power spectum array
* \param fft_size_2 is half of the fft size
* \param masking_thresholds the masking thresholds obtained in db scale
* \param spreaded_unity_gain_bark_spectrum correction to be applied to SSF convolution
* \param spl_reference_values defines the reference values for each bin to convert from db to db SPL
*/
static void
compute_masking_thresholds(float *bark_z, float *absolute_thresholds, float *SSF,
                           float *spectrum, int fft_size_2, float *masking_thresholds,
                           float *spreaded_unity_gain_bark_spectrum,
                           float *spl_reference_values)
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
    tonality_factor = compute_tonality_factor(spectrum, intermediate_band_bins, n_bins_per_band, j); //Uses power spectrum

    //Masking offset
    masking_offset[j] = (tonality_factor * (14.5f + (float)(j + 1)) + 5.5f * (1.f - tonality_factor));

#if BIAS
    //Using offset proposed by Virag (an optimization not needed)
    masking_offset[j] = relative_thresholds[j];
    //Consider tonal noise in upper bands (j>15) due to musical noise of the power Sustraction that was used at First
    if (j > 15)
      masking_offset[j] += HIGH_FREQ_BIAS;
#endif

    //spread Masking threshold
    threshold_j[j] = powf(10.f, log10f(spreaded_spectrum[j]) - (masking_offset[j] / 10.f));

    //Renormalization
    threshold_j[j] -= 10.f * log10f(spreaded_unity_gain_bark_spectrum[j]);

    //Relating the spread masking threshold to the critical band masking thresholds
    //Border case
    if (j == 0)
    {
      start_pos = 0;
    }
    else
    {
      start_pos = intermediate_band_bins[j - 1];
    }
    end_pos = intermediate_band_bins[j];

    for (k = start_pos; k < end_pos; k++)
    {
      masking_thresholds[k] = threshold_j[j];
    }
  }

  //Masking thresholds need to be converted to db spl scale in order to be compared with
  //absolute threshold of hearing
  convert_to_dbspl(spl_reference_values, masking_thresholds, fft_size_2);

  //Take into account the absolute_thresholds of hearing
  for (k = 0; k <= fft_size_2; k++)
  {
    masking_thresholds[k] = MAX(masking_thresholds[k], absolute_thresholds[k]);
  }
}

/**
* alpha and beta computation according to Virags paper. Alphas refers to the oversubtraction
* factor for each fft bin and beta to the spectral flooring. What the oversubtraction
* factor for each bin really does is scaling the noise profile in order reduce more or
* less noise in the supression rule. Spectral flooring limits the amount of reduction in
* each bin. Using masking thresholds for means of adapting this two parameters correlate
* much more with human hearing and results in smoother results than using non linear
* subtraction or others methods of adapting them. Spectral flooring is not used since
* users decide the amount of noise reduccion themselves and spectral flooring is tied to
* that parameter instead of being setted automatically.
* \param fft_p2 the power spectrum of current frame
* \param noise_thresholds_p2 the noise thresholds for each bin estimated previously
* \param fft_size_2 is half of the fft size
* \param alpha_masking is the array of oversubtraction factors for each bin
* \param beta_masking is the array of the spectral flooring factors for each bin
* \param bark_z defines the bark to linear mapping for current spectrum config
* \param absolute_thresholds defines the absolute thresholds of hearing for current spectrum config
* \param SSF defines the spreading function matrix
* \param spreaded_unity_gain_bark_spectrum correction to be applied to SSF convolution
* \param spl_reference_values defines the reference values for each bin to convert from db to db SPL
* \param masking_value is the limit max oversubtraction to be computed
* \param reduction_value is the limit max the spectral flooring to be computed
*/
static void
compute_alpha_and_beta(float *fft_p2, float *noise_thresholds_p2, int fft_size_2,
                       float *alpha_masking, float *beta_masking, float *bark_z,
                       float *absolute_thresholds, float *SSF,
                       float *spreaded_unity_gain_bark_spectrum,
                       float *spl_reference_values, float masking_value,
                       float reduction_value)
{
  int k;
  float masking_thresholds[fft_size_2 + 1];
  float estimated_clean[fft_size_2 + 1];
  float normalized_value;

  //Noise masking threshold must be computed from a clean signal
  //therefor we aproximate a clean signal using a power Sustraction over
  //the original noisy one

  //basic Power Sustraction to estimate clean signal
  for (k = 0; k <= fft_size_2; k++)
  {
    estimated_clean[k] = MAX(fft_p2[k] - noise_thresholds_p2[k], FLT_MIN);
  }

  //Now we can compute noise masking threshold from this clean signal
  compute_masking_thresholds(bark_z, absolute_thresholds, SSF, estimated_clean,
                             fft_size_2, masking_thresholds,
                             spreaded_unity_gain_bark_spectrum, spl_reference_values);

  //First we need the maximun and the minimun value of the masking threshold
  float max_masked_tmp = max_spectral_value(masking_thresholds, fft_size_2);
  float min_masked_tmp = min_spectral_value(masking_thresholds, fft_size_2);

  for (k = 0; k <= fft_size_2; k++)
  {
    //new alpha and beta vector
    if (masking_thresholds[k] == max_masked_tmp)
    {
      alpha_masking[k] = ALPHA_MIN;
      beta_masking[k] = BETA_MIN;
    }
    if (masking_thresholds[k] == min_masked_tmp)
    {
      alpha_masking[k] = masking_value;
      beta_masking[k] = reduction_value;
    }
    if (masking_thresholds[k] < max_masked_tmp && masking_thresholds[k] > min_masked_tmp)
    {
      //Linear interpolation of the value between max and min masked threshold values
      normalized_value = (masking_thresholds[k] - min_masked_tmp) / (max_masked_tmp - min_masked_tmp);

      alpha_masking[k] = (1.f - normalized_value) * ALPHA_MIN + normalized_value * masking_value;
      beta_masking[k] = (1.f - normalized_value) * BETA_MIN + normalized_value * reduction_value;
    }
  }
}
