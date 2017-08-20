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

#include <math.h>
#include <float.h>
#include <stdbool.h>

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

//Window types
#define HANN_WINDOW 0
#define HAMMING_WINDOW 1
#define BLACKMAN_WINDOW 2

#define M_PI 3.14159265358979323846f

#define SP_MAX_NUM 100 //Max number of spectral peaks to find
#define SP_THRESH 0.1f //Threshold to discriminate peaks (high value to discard noise) Linear 0<>1
#define SP_USE_P_INTER true //Use parabolic interpolation
#define SP_MAX_FREQ 16000.f //Highest frequency to search for peaks
#define SP_MIN_FREQ 40.f //Lowest frequency to search for peaks

#define SE_RESOLUTION 100.f //Spectral envelope resolution

//struct for spectral peaks array
typedef struct
{
	float magnitude;
	int position;
} FFTPeak;

//Force already-denormal float value to zero
float
sanitize_denormal(float value)
{
  if (isnan(value))
  {
    return FLT_MIN; //to avoid log errors
    //return 0.f; //to avoid log errors
  }
  else
  {
    return value;
  }
}

int
sign(float x)
{
  return (x >= 0.f ? 1.f : -1.f);
}

int
next_pow_two(int x)
{
	int power = 2;
	while(x >>= 1)
		power <<= 1;
	return power;
}

int
nearest_odd(int x)
{
	if(x%2 == 0)
		return x+1;
	else
		return x;
}

int
nearest_even(int x)
{
	if(x%2 == 0)
		return x;
	else
		return x-1;
}

void
initialize_array(float* array, float value,int size)
{
	for(int k=0; k<size;k++)
	{
		array[k] = value;
	}
}

//Parabolic interpolation as explained in  https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html
void
parabolic_interpolation(float left_val, float middle_val, float right_val,
                        int current_bin, float* result_val, int* result_bin)
{
  float delta_x = 0.5 * ((left_val - right_val) / (left_val - 2.f*middle_val + right_val));
  *result_bin = current_bin + (int)delta_x;
  *result_val = middle_val - 0.25 * (left_val - right_val) * delta_x;
}

//-----------dB SCALE-----------

//power scales
float
from_dB(float gdb)
{
  return (expf(gdb/10.f*logf(10.f)));
}

float
to_dB(float g)
{
  return (10.f*log10f(g));
}

//-----------FREQ <> INDEX OR BIN------------

float
bin_to_freq(int i, float samp_rate, int N)
{
  return (float) i * (samp_rate / N / 2.f);
}

int
freq_to_bin(float freq, float samp_rate, int N)
{
  return (int) (freq / (samp_rate / N / 2.f));
}

//---------SPECTRAL OPERATIONS-------------

//verifies if the spectrum is full of zeros
bool
is_empty(float* spectrum, int N)
{
  int k;
  for(k = 0;k <= N; k++)
  {
    if(spectrum[k] > FLT_MIN)
    {
      return false;
    }
  }
  return true;
}

//finds the max value of the spectrum
float
max_spectral_value(float* spectrum, int N)
{
  int k;
  float max = spectrum[0];
  for(k = 0; k <= N; k++)
  {
    max = MAX(spectrum[k],max);
  }
  return max;
}

//finds the min value of the spectrum
float
min_spectral_value(float* spectrum, int N)
{
  int k;
  float min = spectrum[0];
  for(k = 0; k <= N; k++)
  {
    min = MIN(spectrum[k],min);
  }
  return min;
}

//Mean value of a spectrum
float
spectral_mean(float* a,int m)
{
    float sum=0.f;
    for(int i=0; i<=m; i++)
        sum+=a[i];
    return(sum/(float)(m+1));
}

//Sum of all values of a spectrum
float
spectral_addition(float* a,int m)
{
    float sum=0.f;
    for(int i=0; i<=m; i++)
        sum+=a[i];
    return sum;
}

//Median value of a spectrum
float
spectral_median(float* x,int n)
{
    float temp;
    int i, j;
    float tmp[n+1];
    memcpy(tmp,x,sizeof(float)*(n+1));
    // the following two loops sort the array x in ascending order
    for(i=0; i<n; i++)
    {
      for(j=i+1; j<=n; j++)
      {
        if(tmp[j] < tmp[i])
        {
          // swap elements
          temp = tmp[i];
          tmp[i] = tmp[j];
          tmp[j] = temp;
        }
      }
    }

    if(n%2==0)
    {
      // if there is an even number of elements, return mean of the two elements in the middle
      return((tmp[n/2] + tmp[n/2 - 1]) / 2.f);
    }
    else
    {
      // else return the element in the middle
      return tmp[n/2];
    }
}

float
spectral_moda(float* x,int n)
{
  float temp[n];
  int i,j,pos_max;
  float max;

  for(i = 0;i<n; i++)
  {
    temp[i]=0.f;
  }

  for(i=0; i<n; i++)
  {
    for(j=i; j<n; j++)
    {
      if(x[j] == x[i]) temp[i]++;
    }
  }

  max=temp[0];
  pos_max = 0;
  for(i=0; i<n; i++)
  {
    if(temp[i] > max)
    {
      pos_max = i;
      max=temp[i];
    }
  }
  return x[pos_max];
}

void
get_normalized_spectum(float* spectrum,int N)
{
	int k;
	float max_value = max_spectral_value(spectrum,N);
	float min_value = min_spectral_value(spectrum,N);

	//Normalizing the noise print
	for(k = 0 ; k <= N ; k++)
  {
		spectrum[k] = (spectrum[k]-min_value)/(max_value-min_value);
	}
}

float
spectral_flux(float* spectrum,float* spectrum_prev,float N)
{
  int i;
  float spectral_flux = 0.f;
  float temp;

  for(i = 0;i <= N; i++)
  {
    temp = sqrtf(spectrum[i]) - sqrtf(spectrum_prev[i]); //Recieves power spectrum uses magnitude
		if(temp > 0.f)
			spectral_flux += (temp + fabs(temp))/2.f;
  }
  return spectral_flux;
}

float
high_frequency_content(float* spectrum,float* spectrum_prev,float N)
{
  int i;
  float sum = 0.f;

  for(i = 0;i <= N; i++)
  {
    sum += i*spectrum[i];
  }
  return sum/(float)(N+1);
}

void
spectral_whitening(float* spectrum,float b,int N)
{
	// float peaks_magnitude[peaks_count];
	// float peaks_frequencies[peak_count];
	//
	// //Convert input linear magnitudes to dB scale
	// for (k = 0; k < peaks_count; k++) {
	// 	peaks_magnitude[k] = 2.f*to_dB(spectral_peaks->magnitudes[k]);
	// }
	//
	// //get max peak
	// float max_value = max_spectral_value(peaks_magnitude, peaks_count);

  for (int k = 0; k <= N; k++)
  {
    if(spectrum[k] > FLT_MIN)
    {
      spectrum[k] = (1.f - b)*spectrum[k] + b*(1.f - spectrum[k]);
    }
  }
}

void
spectral_envelope(int fft_size_2, float* fft_p2, int samp_rate, float* spectral_envelope_values)
{
	int k;

	//compute envelope
	int spec_size = fft_size_2+1;
	float spectral_range = bin_to_freq(spec_size,samp_rate,fft_size_2*2);
	int hop = (int)freq_to_bin(SE_RESOLUTION,samp_rate,fft_size_2*2);//Experimental

	for (k = 0; k <= fft_size_2; k+=hop)
	{
		float freq = bin_to_freq(k,samp_rate,fft_size_2*2);

		float bf = freq - MAX(50.0, freq * 0.34); // 0.66
		float ef = freq + MAX(50.0, freq * 0.58); // 1.58
		int b = (int)(bf / spectral_range * (spec_size - 1.0) + 0.5);
		int e = (int)(ef / spectral_range * (spec_size - 1.0) + 0.5);
		b = MAX(b, 0);
		b = MIN(spec_size - 1, b);
		e = MAX(e, b + 1);
		e = MIN(spec_size, e);
		float c = b/2.0 + e/2.0;
		float half_window_length = e - c;

		float n = 0.0;
		float wavg = 0.0;

		for (int i = b; i < e; ++i)
		{
			float weight = 1.0 - fabs((float)(i)-c) / half_window_length;
			weight *= weight;
			weight *= weight;
			float spectrum_energy_val = fft_p2[i];// * fft_p2[i];
			weight *= spectrum_energy_val;
			wavg += spectrum_energy_val * weight;
			n += weight;
		}
		if (n != 0.0)
			wavg /= n;

		//final value
		spectral_envelope_values[k] = wavg;//sqrtf(wavg);
	}
}

void
spectral_peaks(int fft_size_2, float* fft_p2, FFTPeak* spectral_peaks, int* peak_pos,
               int* peaks_count, int samp_rate)
{
  int k;
  float fft_magnitude_db[fft_size_2+1];
  float peak_threshold_db = to_dB(SP_THRESH);
  int max_bin = MIN(freq_to_bin(SP_MAX_FREQ,samp_rate,fft_size_2*2),fft_size_2+1);
  int min_bin = MAX(freq_to_bin(SP_MIN_FREQ,samp_rate,fft_size_2*2),0);
  int result_bin;
  float result_val;

  //Get the magnitude spectrum in dB scale (twise as precise than using linear scale)
  for(k = 0; k<=fft_size_2;k++)
  {
    fft_magnitude_db[k] = to_dB(sqrtf(fft_p2[k]));
  }

  //index for the magnitude array
  int i = min_bin;

  //Index for peak array
  k = 0;

  //The zero bin could be a peak
  if (i+1 < fft_size_2+1 && fft_magnitude_db[i] > fft_magnitude_db[i+1])
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
  while(k < SP_MAX_NUM || i < max_bin)
  {
    //descending a peak
    while (i+1 < fft_size_2 && fft_magnitude_db[i] >= fft_magnitude_db[i+1])
    {
      i++;
    }
    //ascending a peak
    while (i+1 < fft_size_2 && fft_magnitude_db[i] < fft_magnitude_db[i+1])
    {
      i++;
    }

    //when reaching a peak verify that is one value peak or multiple values peak
    int j = i;
    while (j+1 < fft_size_2 && (fft_magnitude_db[j] == fft_magnitude_db[j+1]))
    {
      j++;
    }

    //end of the flat peak if the peak decreases is really a peak otherwise it is not
    if (j+1 < fft_size_2 && fft_magnitude_db[j+1] < fft_magnitude_db[j] && fft_magnitude_db[j] > peak_threshold_db)
    {
      result_bin = 0.0;
      result_val = 0.0;

      if (j != i) { //peak between i and j
        if (SP_USE_P_INTER)
        {
          result_bin = (i + j) * 0.5;//center bin of the flat peak
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
          parabolic_interpolation(fft_magnitude_db[j-1], fft_magnitude_db[j], fft_magnitude_db[j+1], j, &result_val, &result_bin);
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
    if (i+1 >= fft_size_2)
    {
      if (i == fft_size_2-1 && fft_magnitude_db[i-1] < fft_magnitude_db[i] &&
        fft_magnitude_db[i+1] < fft_magnitude_db[i] &&
        fft_magnitude_db[i] > peak_threshold_db)
        {
          result_bin = 0.0;
          result_val = 0.0;
          if (SP_USE_P_INTER)
          {
            parabolic_interpolation(fft_magnitude_db[i-1], fft_magnitude_db[i], fft_magnitude_db[i+1], j, &result_val, &result_bin);
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

//---------------TIME SMOOTHING--------------

//This was proposed in this work SPECTRAL SUBTRACTION WITH ADAPTIVE AVERAGING OF THE GAIN FUNCTION
void
spectrum_adaptive_time_smoothing(int fft_size_2, float* spectrum_prev, float* spectrum,
                                 float* noise_thresholds, float* prev_beta, float coeff)
{
  int k;
  float discrepancy, numerator = 0.f, denominator = 0.f;
  float beta_ts;
  float beta_smooth;
  float gamma_ts;

  for (k = 0; k <= fft_size_2; k++)
  {
    //These has to be magnitude spectrums
    numerator += fabs(spectrum[k] - noise_thresholds[k]);
    denominator += noise_thresholds[k];
  }
  //this is the discrepancy of the spectum
  discrepancy = numerator/denominator;
  //beta is the adaptive coefficient
  beta_ts = MIN(discrepancy,1.f);

  //printf("%f\n", beta_ts);

  //Gamma is the smoothing coefficient of the adaptive factor beta
  if(*prev_beta < beta_ts)
  {
    gamma_ts = 0.f;
  }
  else
  {
    gamma_ts = coeff;
  }

  //Smoothing beta
  beta_smooth = gamma_ts * *(prev_beta) + (1.f - gamma_ts)*beta_ts;

  //copy current value to previous
  *prev_beta = beta_smooth;

  //Apply the adaptive smoothed beta over the signal
  for (k = 0; k <= fft_size_2; k++)
  {
    spectrum[k] = (1.f - beta_smooth) * spectrum_prev[k] + beta_smooth * spectrum[k];
  }
}

void
apply_time_envelope(float* spectrum, float* spectrum_prev, float N, float release_coeff)
{
  int k;

  for (k = 0; k <= N ; k++)
  {
    //It doesn't make much sense to have an attack slider when there is time smoothing
    if (spectrum[k] > spectrum_prev[k])
    {
      //Release (when signal is incrementing in amplitude)
      spectrum[k] = release_coeff*spectrum_prev[k] + (1.f-release_coeff)*spectrum[k];
    }
	}
}

//-----------WINDOW---------------

//blackman window values computing
float
blackman(int k, int N)
{
  float p = ((float)(k))/((float)(N));
  return 0.42-0.5*cosf(2.f*M_PI*p) + 0.08*cosf(4.f*M_PI*p);
}

//hanning window values computing
float
hanning(int k, int N)
{
  float p = ((float)(k))/((float)(N));
  return 0.5 - 0.5 * cosf(2.f*M_PI*p);
}

//hamming window values computing
float
hamming(int k, int N)
{
  float p = ((float)(k))/((float)(N));
  return 0.54 - 0.46 * cosf(2.f*M_PI*p);
}

//wrapper to compute windows values
void
fft_window(float* window, int N, int window_type)
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
    }
  }
}

float
get_window_scale_factor(float* input_window,float* output_window,int frame_size)
{
 	float sum=0.f;
 	for(int i=0; i<frame_size; i++)
 			sum+= input_window[i]*output_window[i];
 	return (sum/(float)(frame_size));
}

//wrapper for pre and post processing windows
void
fft_pre_and_post_window(float* input_window, float* output_window, int frame_size,
                        int window_option_input, int window_option_output,
                        float* overlap_scale_factor)
{
  //Input window
  switch(window_option_input)
  {
    case 0: // HANN-HANN
      fft_window(input_window,frame_size,0); //STFT input window
      break;
    case 1: //HAMMING-HANN
      fft_window(input_window,frame_size,1); //STFT input window
      break;
    case 2: //BLACKMAN-HANN
      fft_window(input_window,frame_size,2); //STFT input window
      break;
  }

  //Output window
  switch(window_option_output)
  {
    case 0: // HANN-HANN
      fft_window(output_window,frame_size,0); //STFT input window
      break;
    case 1: //HAMMING-HANN
      fft_window(output_window,frame_size,1); //STFT input window
      break;
    case 2: //BLACKMAN-HANN
      fft_window(output_window,frame_size,2); //STFT input window
      break;
  }

  //Scaling necessary for perfect reconstruction using Overlapp Add
  *(overlap_scale_factor) = get_window_scale_factor(input_window,output_window,frame_size);
}

void
get_info_from_bins(float* fft_p2, float* fft_magnitude, float* fft_phase,
									 int fft_size_2, int fft_size, float* fft_buffer)
{
	int k;
	float real_p,imag_n,mag,p2,phase;

	//Look at http://www.fftw.org/fftw2_doc/fftw_2.html
	//DC bin
	real_p = fft_buffer[0];
	imag_n = 0.f;

	fft_p2[0] = real_p*real_p;
	fft_magnitude[0] = real_p;
	fft_phase[0] = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist

	//Get the rest of positive spectrum and compute the magnitude
	for (k = 1; k <= fft_size_2; k++)
	{
		//Get the half complex spectrum reals and complex
		real_p = fft_buffer[k];
		imag_n = fft_buffer[fft_size - k];

		//Get the magnitude, phase and power spectrum
		if(k < fft_size_2)
		{
			p2 = (real_p*real_p + imag_n*imag_n);
			mag = sqrtf(p2);//sqrt(real^2+imag^2)
			phase = atan2f(real_p, imag_n);
		}
		else
		{
			//Nyquist - this is due to half complex transform look at http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html
			p2 = real_p*real_p;
			mag = real_p;
			phase = atan2f(real_p, 0.f); //Phase is 0 for DC and nyquist
		}
		//Store values in magnitude and power arrays (this stores the positive spectrum only)
		fft_p2[k] = p2;
		fft_magnitude[k] = mag; //This is not used but part of the STFT transform for generic use
		fft_phase[k] = phase; //This is not used but part of the STFT transform for generic use
	}
}
