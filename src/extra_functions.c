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
* \file extra_functions.c
* \author Luciano Dato
* \brief Extra methods used by others. This keeps clean other files.
*/

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

//Window types
#define HANN_WINDOW 0
#define HAMMING_WINDOW 1
#define BLACKMAN_WINDOW 2
#define VORBIS_WINDOW 3

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#define SP_MAX_NUM 100      //Max number of spectral peaks to find
#define SP_THRESH 0.1f      //Threshold to discriminate peaks (high value to discard noise) Linear 0<>1
#define SP_USE_P_INTER true //Use parabolic interpolation
#define SP_MAX_FREQ 16000.f //Highest frequency to search for peaks
#define SP_MIN_FREQ 40.f    //Lowest frequency to search for peaks

#define SE_RESOLUTION 100.f //Spectral envelope resolution

#define TP_UPPER_LIMIT 5.f //This correspond to the upper limit of the adaptive threshold multiplier. Should be the same as the ttl configured one

/**
* Method to force already-denormal float value to zero.
* \param value to sanitize
*/
float sanitize_denormal(float value)
{
  if (isnan(value))
  {
    return FLT_MIN; //to avoid log errors
  }
  else
  {
    return value;
  }
}

///sign function.
int sign(float x)
{
  return (x >= 0.f ? 1.f : -1.f);
}

///gets the next power of two of a number x.
int next_pow_two(int x)
{
  int power = 2;
  while (x >>= 1)
    power <<= 1;
  return power;
}

///gets the nearest odd number of a number x.
int nearest_odd(int x)
{
  if (x % 2 == 0)
    return x + 1;
  else
    return x;
}

///gets the nearest even number of a number x.
int nearest_even(int x)
{
  if (x % 2 == 0)
    return x;
  else
    return x - 1;
}

///converts a db value to linear scale.
float from_dB(float gdb)
{
  return (expf(gdb / 10.f * logf(10.f)));
}

///converts a linear value to db scale.
float to_dB(float g)
{
  return (10.f * log10f(g));
}

/*Maps a bin number to a frequency
* \param i bin number
* \param samp_rate current sample rate of the host
* \param N size of the fft
*/
float bin_to_freq(int i, float samp_rate, int N)
{
  return (float)i * (samp_rate / N / 2.f);
}

/*Maps a frequency to a bin number
* \param freq frequency
* \param samp_rate current sample rate of the host
* \param N size of the fft
*/
int freq_to_bin(float freq, float samp_rate, int N)
{
  return (int)(freq / (samp_rate / N / 2.f));
}

//---------SPECTRAL METHODS-------------

/**
* FFT peak struct. To describe spectral peaks.
*/
typedef struct
{
  float magnitude;
  int position;
} FFTPeak;

/**
* Parabolic interpolation as explained in  https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html.
* This is used for more precise spectral peak detection.
* \param left_val value at the left of the point to interpolate
* \param middle_val value at the middle of the point to interpolate
* \param right_val value at the right of the point to interpolate
* \param current_bin current bin value before interpolation
* \param result_val interpolation value result
* \param result_val interpolation bin result
*/
void parabolic_interpolation(float left_val, float middle_val, float right_val,
                             int current_bin, float *result_val, int *result_bin)
{
  float delta_x = 0.5 * ((left_val - right_val) / (left_val - 2.f * middle_val + right_val));
  *result_bin = current_bin + (int)delta_x;
  *result_val = middle_val - 0.25 * (left_val - right_val) * delta_x;
}

/**
* To initialize an array to a single value in all positions.
* \param array the array to initialize
* \param value the value to copy to every position in the array
* \param size the size of the array
*/
void initialize_array(float *array, float value, int size)
{
  for (int k = 0; k < size; k++)
  {
    array[k] = value;
  }
}

/**
* Verifies if the spectrum is full of zeros.
* \param spectrum the array to check
* \param N the size of the array (half the fft size plus 1)
*/
bool is_empty(float *spectrum, int N)
{
  int k;
  for (k = 0; k <= N; k++)
  {
    if (spectrum[k] > FLT_MIN)
    {
      return false;
    }
  }
  return true;
}

/**
* Finds the max value of the spectrum.
* \param spectrum the array to check
* \param N the size of the array (half the fft size plus 1)
*/
float max_spectral_value(float *spectrum, int N)
{
  int k;
  float max = spectrum[0];
  for (k = 0; k <= N; k++)
  {
    max = MAX(spectrum[k], max);
  }
  return max;
}

/**
* Finds the min value of the spectrum.
* \param spectrum the array to check
* \param N the size of the array (half the fft size plus 1)
*/
float min_spectral_value(float *spectrum, int N)
{
  int k;
  float min = spectrum[0];
  for (k = 0; k <= N; k++)
  {
    min = MIN(spectrum[k], min);
  }
  return min;
}

/**
* Finds the mean value of the spectrum.
* \param a the array to check
* \param m the size of the array (half the fft size plus 1)
*/
float spectral_mean(float *a, int m)
{
  float sum = 0.f;
  for (int i = 0; i <= m; i++)
    sum += a[i];
  return (sum / (float)(m + 1));
}

/**
* Sums of all values of a spectrum.
* \param a the array to sum
* \param m the size of the array (half the fft size plus 1)
*/
float spectral_addition(float *a, int m)
{
  float sum = 0.f;
  for (int i = 0; i <= m; i++)
    sum += a[i];
  return sum;
}

/**
* Finds the median value of the spectrum.
* \param x the array to check
* \param n the size of the array (half the fft size plus 1)
*/
float spectral_median(float *x, int n)
{
  float temp;
  int i, j;
  float tmp[n + 1];
  memcpy(tmp, x, sizeof(float) * (n + 1));
  // the following two loops sort the array x in ascending order
  for (i = 0; i < n; i++)
  {
    for (j = i + 1; j <= n; j++)
    {
      if (tmp[j] < tmp[i])
      {
        // swap elements
        temp = tmp[i];
        tmp[i] = tmp[j];
        tmp[j] = temp;
      }
    }
  }

  if (n % 2 == 0)
  {
    // if there is an even number of elements, return mean of the two elements in the middle
    return ((tmp[n / 2] + tmp[n / 2 - 1]) / 2.f);
  }
  else
  {
    // else return the element in the middle
    return tmp[n / 2];
  }
}

/**
* Finds the moda value of the spectrum.
* \param x the array to check
* \param n the size of the array (half the fft size plus 1)
*/
float spectral_moda(float *x, int n)
{
  float temp[n];
  int i, j, pos_max;
  float max;

  for (i = 0; i < n; i++)
  {
    temp[i] = 0.f;
  }

  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      if (x[j] == x[i])
        temp[i]++;
    }
  }

  max = temp[0];
  pos_max = 0;
  for (i = 0; i < n; i++)
  {
    if (temp[i] > max)
    {
      pos_max = i;
      max = temp[i];
    }
  }
  return x[pos_max];
}

/**
* Normalizes spectral values.
* \param spectrum the spectrum to normalize
* \param N the size of the spectrum (half the fft size plus 1)
*/
void get_normalized_spectum(float *spectrum, int N)
{
  int k;
  float max_value = max_spectral_value(spectrum, N);
  float min_value = min_spectral_value(spectrum, N);

  //Normalizing the noise print
  for (k = 0; k <= N; k++)
  {
    spectrum[k] = (spectrum[k] - min_value) / (max_value - min_value);
  }
}

/**
* Outputs the spectral flux between two spectrums.
* \param spectrum the current power spectrum
* \param spectrum_prev the previous power spectrum
* \param N the size of the spectrum (half the fft size plus 1)
*/
float spectral_flux(float *spectrum, float *spectrum_prev, float N)
{
  int i;
  float spectral_flux = 0.f;
  float temp;

  for (i = 0; i <= N; i++)
  {
    temp = sqrtf(spectrum[i]) - sqrtf(spectrum_prev[i]); //Recieves power spectrum uses magnitude
    spectral_flux += (temp + fabs(temp)) / 2.f;
  }
  return spectral_flux;
}

/**
* Outputs the high frequency content of the spectrum.
* \param spectrum the current power spectrum
* \param N the size of the spectrum (half the fft size plus 1)
*/
float high_frequency_content(float *spectrum, float N)
{
  int i;
  float sum = 0.f;

  for (i = 0; i <= N; i++)
  {
    sum += i * spectrum[i];
  }
  return sum / (float)(N + 1);
}

/**
* Computes the spectral envelope like Robel 'Efficient Spectral Envelope Estimation and its
* application to pitch shifting and envelope preservation' indicates.
* \param fft_size_2 half of the fft size
* \param fft_p2 the current power spectrum
* \param samp_rate current sample rate of the host
* \param spectral_envelope_values array that holds the spectral envelope values
*/
void spectral_envelope(int fft_size_2, float *fft_p2, int samp_rate, float *spectral_envelope_values)
{
  int k;

  //compute envelope
  int spec_size = fft_size_2 + 1;
  float spectral_range = bin_to_freq(spec_size, samp_rate, fft_size_2 * 2);
  int hop = (int)freq_to_bin(SE_RESOLUTION, samp_rate, fft_size_2 * 2); //Experimental

  for (k = 0; k <= fft_size_2; k += hop)
  {
    float freq = bin_to_freq(k, samp_rate, fft_size_2 * 2);

    float bf = freq - MAX(50.0, freq * 0.34); // 0.66
    float ef = freq + MAX(50.0, freq * 0.58); // 1.58
    int b = (int)(bf / spectral_range * (spec_size - 1.0) + 0.5);
    int e = (int)(ef / spectral_range * (spec_size - 1.0) + 0.5);
    b = MAX(b, 0);
    b = MIN(spec_size - 1, b);
    e = MAX(e, b + 1);
    e = MIN(spec_size, e);
    float c = b / 2.0 + e / 2.0;
    float half_window_length = e - c;

    float n = 0.0;
    float wavg = 0.0;

    for (int i = b; i < e; ++i)
    {
      float weight = 1.0 - fabs((float)(i)-c) / half_window_length;
      weight *= weight;
      weight *= weight;
      float spectrum_energy_val = fft_p2[i]; // * fft_p2[i];
      weight *= spectrum_energy_val;
      wavg += spectrum_energy_val * weight;
      n += weight;
    }
    if (n != 0.0)
      wavg /= n;

    //final value
    spectral_envelope_values[k] = wavg; //sqrtf(wavg);
  }
}

/**
* Finds the spectral peaks of a spectrum. (Not used in current version)
* \param fft_size_2 half of the fft size
* \param fft_p2 the current power spectrum
* \param spectral_peaks array of resulting spectral peaks founded
* \param peak_pos array of positions of resulting spectral peaks
* \param peaks_count counter of peaks founded
* \param samp_rate current sample rate of the host
*/
void spectral_peaks(int fft_size_2, float *fft_p2, FFTPeak *spectral_peaks, int *peak_pos,
                    int *peaks_count, int samp_rate)
{
  int k;
  float fft_magnitude_db[fft_size_2 + 1];
  float peak_threshold_db = to_dB(SP_THRESH);
  int max_bin = MIN(freq_to_bin(SP_MAX_FREQ, samp_rate, fft_size_2 * 2), fft_size_2 + 1);
  int min_bin = MAX(freq_to_bin(SP_MIN_FREQ, samp_rate, fft_size_2 * 2), 0);
  int result_bin;
  float result_val;

  //Get the magnitude spectrum in dB scale (twise as precise than using linear scale)
  for (k = 0; k <= fft_size_2; k++)
  {
    fft_magnitude_db[k] = to_dB(sqrtf(fft_p2[k]));
  }

  //index for the magnitude array
  int i = min_bin;

  //Index for peak array
  k = 0;

  //The zero bin could be a peak
  if (i + 1 < fft_size_2 + 1 && fft_magnitude_db[i] > fft_magnitude_db[i + 1])
  {
    if (fft_magnitude_db[i] > peak_threshold_db)
    {
      spectral_peaks[k].position = i;
      spectral_peaks[k].magnitude = sqrtf(from_dB(fft_magnitude_db[i]));
      peak_pos[i] = 1;
      k++;
    }
  }

  //Peak finding loop
  while (k < SP_MAX_NUM || i < max_bin)
  {
    //descending a peak
    while (i + 1 < fft_size_2 && fft_magnitude_db[i] >= fft_magnitude_db[i + 1])
    {
      i++;
    }
    //ascending a peak
    while (i + 1 < fft_size_2 && fft_magnitude_db[i] < fft_magnitude_db[i + 1])
    {
      i++;
    }

    //when reaching a peak verify that is one value peak or multiple values peak
    int j = i;
    while (j + 1 < fft_size_2 && (fft_magnitude_db[j] == fft_magnitude_db[j + 1]))
    {
      j++;
    }

    //end of the flat peak if the peak decreases is really a peak otherwise it is not
    if (j + 1 < fft_size_2 && fft_magnitude_db[j + 1] < fft_magnitude_db[j] && fft_magnitude_db[j] > peak_threshold_db)
    {
      result_bin = 0.0;
      result_val = 0.0;

      if (j != i)
      { //peak between i and j
        if (SP_USE_P_INTER)
        {
          result_bin = (i + j) * 0.5; //center bin of the flat peak
        }
        else
        {
          result_bin = i;
        }
        result_val = fft_magnitude_db[i];
      }
      else
      { //interpolate peak at i-1, i and i+1
        if (SP_USE_P_INTER)
        {
          parabolic_interpolation(fft_magnitude_db[j - 1], fft_magnitude_db[j], fft_magnitude_db[j + 1], j, &result_val, &result_bin);
        }
        else
        {
          result_bin = j;
          result_val = fft_magnitude_db[j];
        }
      }

      spectral_peaks[k].position = result_bin;
      spectral_peaks[k].magnitude = sqrtf(from_dB(result_val));
      peak_pos[i] = 1;
      k++;
    }

    //if turned out not to be a peak advance i
    i = j;

    //If it's the last position of the array
    if (i + 1 >= fft_size_2)
    {
      if (i == fft_size_2 - 1 && fft_magnitude_db[i - 1] < fft_magnitude_db[i] &&
          fft_magnitude_db[i + 1] < fft_magnitude_db[i] &&
          fft_magnitude_db[i] > peak_threshold_db)
      {
        result_bin = 0.0;
        result_val = 0.0;
        if (SP_USE_P_INTER)
        {
          parabolic_interpolation(fft_magnitude_db[i - 1], fft_magnitude_db[i], fft_magnitude_db[i + 1], j, &result_val, &result_bin);
        }
        else
        {
          result_bin = i;
          result_val = fft_magnitude_db[i];
        }
        spectral_peaks[k].position = result_bin;
        spectral_peaks[k].magnitude = sqrtf(from_dB(result_val));
        peak_pos[i] = 1;
        k++;
      }
      break;
    }
  }
  *peaks_count = k;
  //printf("%i\n",k );
}

/**
* Outputs the p norm of a spectrum.
* \param spectrum the power spectrum to get the norm of
* \param N the size of the array
* \param p the norm number
*/
float spectrum_p_norm(float *spectrum, float N, float p)
{
  float sum = 0.f;

  for (int k = 0; k < N; k++)
  {
    sum += powf(spectrum[k], p);
  }

  return powf(sum, 1.f / p);
}

//-----------WINDOW---------------

/**
* blackman window values computing.
* \param k bin number
* \param N fft size
*/
float blackman(int k, int N)
{
  float p = ((float)(k)) / ((float)(N));
  return 0.42 - 0.5 * cosf(2.f * M_PI * p) + 0.08 * cosf(4.f * M_PI * p);
}

/**
* hanning window values computing.
* \param k bin number
* \param N fft size
*/
float hanning(int k, int N)
{
  float p = ((float)(k)) / ((float)(N));
  return 0.5 - 0.5 * cosf(2.f * M_PI * p);
}

/**
* hamming window values computing.
* \param k bin number
* \param N fft size
*/
float hamming(int k, int N)
{
  float p = ((float)(k)) / ((float)(N));
  return 0.54 - 0.46 * cosf(2.f * M_PI * p);
}

/**
* Vorbis window values computing. It satisfies Princen-Bradley criterion so perfect
* reconstruction could be achieved with 50% overlap when used both in Analysis and
* Synthesis
* \param k bin number
* \param N fft size
*/
float vorbis(int k, int N)
{
  float p = ((float)(k)) / ((float)(N));
  return sinf(M_PI / 2.f * powf(sinf(M_PI * p), 2.f));
}

/**
* Wrapper to compute windows values.
* \param window array for window values
* \param N fft size
* \param window_type type of window
*/
void fft_window(float *window, int N, int window_type)
{
  int k;
  for (k = 0; k < N; k++)
  {
    switch (window_type)
    {
    case BLACKMAN_WINDOW:
      window[k] = blackman(k, N);
      break;
    case HANN_WINDOW:
      window[k] = hanning(k, N);
      break;
    case HAMMING_WINDOW:
      window[k] = hamming(k, N);
      break;
    case VORBIS_WINDOW:
      window[k] = vorbis(k, N);
      break;
    }
  }
}

/**
* Gets the magnitude and phase spectrum of the complex spectrum. Takimg into account that
* the half complex fft was used half of the spectrum contains the real part the other
* the imaginary. Look at http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html for
* more info. DC bin was treated as suggested in http://www.fftw.org/fftw2_doc/fftw_2.html
* \param fft_p2 the current power spectrum
* \param fft_magnitude the current magnitude spectrum
* \param fft_phase the current phase spectrum
* \param fft_size_2 half of the fft size
* \param fft_size size of the fft
* \param fft_buffer buffer with the complex spectrum of the fft transform
*/
void get_info_from_bins(float *fft_p2, float *fft_magnitude, float *fft_phase,
                        int fft_size_2, int fft_size, float *fft_buffer)
{
  int k;
  float real_p, imag_n, mag, p2, phase;

  //DC bin
  real_p = fft_buffer[0];
  imag_n = 0.f;

  fft_p2[0] = real_p * real_p;
  fft_magnitude[0] = real_p;
  fft_phase[0] = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist

  //Get the rest of positive spectrum and compute the magnitude
  for (k = 1; k <= fft_size_2; k++)
  {
    //Get the half complex spectrum reals and complex
    real_p = fft_buffer[k];
    imag_n = fft_buffer[fft_size - k];

    //Get the magnitude, phase and power spectrum
    if (k < fft_size_2)
    {
      p2 = (real_p * real_p + imag_n * imag_n);
      mag = sqrtf(p2); //sqrt(real^2+imag^2)
      phase = atan2f(real_p, imag_n);
    }
    else
    {
      //Nyquist - this is due to half complex transform
      p2 = real_p * real_p;
      mag = real_p;
      phase = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist
    }
    //Store values in magnitude and power arrays (this stores the positive spectrum only)
    fft_p2[k] = p2;
    fft_magnitude[k] = mag; //This is not used but part of the STFT transform for generic use
    fft_phase[k] = phase;   //This is not used but part of the STFT transform for generic use
  }
}
