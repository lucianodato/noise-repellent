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

//Window types
#define HANN_WINDOW 0
#define HAMMING_WINDOW 1
#define BLACKMAN_WINDOW 2

#define M_PI 3.14159265358979323846f

#define ONSET_THRESH 100.f //For onset detection

#define SP_MAX_NUM 100 //Max number of spectral peaks to find
#define SP_THRESH 0.1f //Threshold to discriminate peaks (high value to discard noise) Linear 0<>1
#define SP_USE_P_INTER 1 //Use parabolic interpolation
#define SP_MAX_FREQ 16000.f //Highest frequency to search for peaks
#define SP_MIN_FREQ 40.f //Lowest frequency to search for peaks

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

//This is from wikipedia ;)
float savgol_quad_5[5] = {-0.085714,0.342857,0.485714,0.342857,-0.085714};
float savgol_quad_7[7] = {-0.095238,0.142857,0.285714,0.333333,0.285714,0.142857,-0.095238};
float savgol_quad_9[9] = {-0.090909,0.060606,0.168831,0.233766,0.255411,0.233766,0.168831,0.060606,-0.090909};
float savgol_quad_11[11] = {-0.083916,0.020979,0.102564,0.160839,0.195804,0.207459,0.195804,0.160839,0.102564,0.020979,-0.083916};
float savgol_quad_13[13] = {-0.076923,0.000000,0.062937,0.111888,0.146853,0.167832,0.174825,0.167832,0.146853,0.111888,0.062937,0.000000,-0.076923};
float savgol_quad_15[15] = {-0.070588,-0.011765,0.038009,0.078733,0.110407,0.133032,0.146606,0.151131,0.146606,0.133032,0.110407,0.078733,0.038009,-0.011765,-0.070588};
float savgol_quad_17[17] = {-0.065015,-0.018576,0.021672,0.055728,0.083591,0.105263,0.120743,0.130031,0.133127,0.130031,0.120743,0.105263,0.083591,0.055728,0.021672,-0.018576,-0.065015};
float savgol_quad_19[19] = {-0.060150,-0.022556,0.010615,0.039363,0.063689,0.083591,0.099071,0.110128,0.116762,0.118974,0.116762,0.110128,0.099071,0.083591,0.063689,0.039363,0.010615,-0.022556,-0.060150};
float savgol_quad_21[21] = {-0.055901,-0.024845,0.002942,0.027460,0.048709,0.066688,0.081399,0.092841,0.101013,0.105917,0.107551,0.105917,0.101013,0.092841,0.081399,0.066688,0.048709,0.027460,0.002942,-0.024845,-0.055901};
float savgol_quad_23[23] = {-0.052174,-0.026087,-0.002484,0.018634,0.037267,0.053416,0.067081,0.078261,0.086957,0.093168,0.096894,0.098137,0.096894,0.093168,0.086957,0.078261,0.067081,0.053416,0.037267,0.018634,-0.002484,-0.026087,-0.052174};
float savgol_quad_25[25] = {-0.048889,-0.026667,-0.006377,0.011981,0.028406,0.042899,0.055459,0.066280,0.074783,0.081546,0.086377,0.089275,0.090242,0.089275,0.086377,0.081546,0.074783,0.066280,0.055459,0.042899,0.028406,0.011981,-0.006377,-0.026667,-0.048889};
float savgol_quart_7[7] = {0.021645,-0.129870,0.324675,0.567100,0.324675,-0.129870,0.021645};
float savgol_quart_9[9] = {0.034965,-0.128205,0.069930,0.314685,0.417249,0.314685,0.069930,-0.128205,0.034965};

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

//---------------WHITENING--------------

void
whitening(float* spectrum,float b,int N)
{
  for (int k = 0; k < N; k++)
  {
    if(spectrum[k] > FLT_MIN)
    {
      spectrum[k] = (1.f - b)*spectrum[k] + b*(1.f - spectrum[k]);
    }
  }
}

//---------------TRANSIENTS--------------

float
spectral_flux(float* spectrum,float* spectrum_prev,float N)
{
  int i;
  float spectral_flux = 0.f;
  float temp;

  for(i = 0;i <= N; i++)
  {
    temp = sqrtf(spectrum[i]) - sqrtf(spectrum_prev[i]); //Recieves power spectrum uses magnitude
    spectral_flux += (temp + fabs(temp))/2.f;
  }
  return spectral_flux;
}

float
transient_preservation(float* spectrum,float* spectrum_prev,float N)
{
  float spectral_flux_value = spectral_flux(spectrum, spectrum_prev, N);

  if (spectral_flux_value > ONSET_THRESH) //This is poor sounding maybe the best approch is multiresolution TODO
    return 1.f/spectral_flux_value;
  else
    return 1.f;
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
    numerator += fabs(sqrtf(spectrum[k]) - sqrtf(noise_thresholds[k]));
    denominator += sqrtf(noise_thresholds[k]);
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
    spectrum[k] = (1.f - beta_ts) * spectrum[k] + beta_ts * spectrum_prev[k];
  }
}

void
apply_envelope(float* spectrum, float* spectrum_prev, float N, float release_coeff)
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

//---------------------SPECTRAL PEAK DETECTION--------------------------

//Peak interpolation based on parabolic curve as explained in https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html
void
parabolic_interpolation(float left_val, float middle_val, float right_val,
                        int current_bin, float* result_val, int* result_bin)
{
  float delta_x = 0.5 * ((left_val - right_val) / (left_val - 2.f*middle_val + right_val));
  *result_bin = current_bin + (int)delta_x;
  *result_val = middle_val - 0.25 * (left_val - right_val) * delta_x;
}

void
spectral_peaks(int fft_size_2, float* fft_p2, FFTPeak* spectral_peaks, int* peak_pos,
               int* peaks_count, int samp_rate)
{
  int k;
  float fft_magnitude_dB[fft_size_2+1];
  float peak_threshold_db = to_dB(SP_THRESH);
  int max_bin = MIN(freq_to_bin(SP_MAX_FREQ,samp_rate,fft_size_2),fft_size_2+1);
  int min_bin = MAX(freq_to_bin(SP_MIN_FREQ,samp_rate,fft_size_2),0);
  int result_bin;
  float result_val;

  //Get the magnitude spectrum in dB scale (twise as precise than using linear scale)
  for(k = 0; k<=fft_size_2;k++)
  {
    fft_magnitude_dB[k] = to_dB(sqrtf(fft_p2[k]));
  }

  //index for the magnitude array
  int i = min_bin;

  //Index for peak array
  k = 0;

  //The zero bin could be a peak
  if (i+1 < fft_size_2+1 && fft_magnitude_dB[i] > fft_magnitude_dB[i+1])
  {
    if (fft_magnitude_dB[i] > peak_threshold_db)
    {
      spectral_peaks[k].position = i;
      spectral_peaks[k].magnitude = sqrtf(from_dB(fft_magnitude_dB[i]));
      peak_pos[i] = 1;
      k++;
    }
  }

  //Peak finding loop
  while(k < SP_MAX_NUM || i < max_bin)
  {
    //descending a peak
    while (i+1 < fft_size_2 && fft_magnitude_dB[i] >= fft_magnitude_dB[i+1])
    {
      i++;
    }
    //ascending a peak
    while (i+1 < fft_size_2 && fft_magnitude_dB[i] < fft_magnitude_dB[i+1])
    {
      i++;
    }

    //when reaching a peak verify that is one value peak or multiple values peak
    int j = i;
    while (j+1 < fft_size_2 && (fft_magnitude_dB[j] == fft_magnitude_dB[j+1]))
    {
      j++;
    }

    //end of the flat peak if the peak decreases is really a peak otherwise it is not
    if (j+1 < fft_size_2 && fft_magnitude_dB[j+1] < fft_magnitude_dB[j] && fft_magnitude_dB[j] > peak_threshold_db)
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
        result_val = fft_magnitude_dB[i];
      }
      else
      { //interpolate peak at i-1, i and i+1
        if (SP_USE_P_INTER)
        {
          parabolic_interpolation(fft_magnitude_dB[j-1], fft_magnitude_dB[j], fft_magnitude_dB[j+1], j, &result_val, &result_bin);
        }
        else
        {
          result_bin = j;
          result_val = fft_magnitude_dB[j];
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
      if (i == fft_size_2-1 && fft_magnitude_dB[i-1] < fft_magnitude_dB[i] &&
        fft_magnitude_dB[i+1] < fft_magnitude_dB[i] &&
        fft_magnitude_dB[i] > peak_threshold_db)
        {
          result_bin = 0.0;
          result_val = 0.0;
          if (SP_USE_P_INTER)
          {
            parabolic_interpolation(fft_magnitude_dB[i-1], fft_magnitude_dB[i], fft_magnitude_dB[i+1], j, &result_val, &result_bin);
          }
          else
          {
            result_bin = i;
            result_val = fft_magnitude_dB[i];
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



//---------------------SPECTRAL SMOOTHERS--------------------------

//Spectral smoothing with rectangular boxcar or unweighted sliding-average smooth
void
spectral_smoothing_MA(float* spectrum, int kernel_width,int N)
{
  int k;
  float smoothing_tmp[N+1];
  float t_spectrum[N+1];

  if (kernel_width == 0) return;

  //Initialize smothingbins_tmp
  for (k = 0; k <= N; ++k)
  {
    t_spectrum[k] = logf(spectrum[k]);
    smoothing_tmp[k] = 0.f;//Initialize temporal spectrum
  }

  for (k = 0; k <= N; ++k)
  {
    const int j0 = MAX(0, k - kernel_width);
    const int j1 = MIN(N, k + kernel_width);
    for(int l = j0; l <= j1; ++l)
    {
      smoothing_tmp[k] += t_spectrum[l];
    }
    smoothing_tmp[k] /= (j1 - j0 + 1);
  }

  for (k = 0; k <= N; ++k)
  {
    spectrum[k] = expf(smoothing_tmp[k]);
    if (k < N)
      spectrum[(2*N)-k] = expf(smoothing_tmp[(2*N)- k]);
  }
}
//Spectral smoothing with median filter
void
spectral_smoothing_MM(float* spectrum, int kernel_width, int N)
{
  int k;
  float smoothing_tmp[N+1];

  if (kernel_width == 0) return;

  for (k = 0; k <= N; ++k)
  {
    const int j0 = MAX(0, k - kernel_width);
    const int j1 = MIN(N, k + kernel_width);

    float aux[j1-j0+1];
    for(int l = j0; l <= j1; ++l)
    {
      aux[l] = spectrum[l];
    }
    smoothing_tmp[k] = spectral_median(aux,j1-j0+1);
  }

  for (k = 0; k <= N; ++k)
  {
    spectrum[k] = smoothing_tmp[k];
  }
}


void
spectral_smoothing_MAH(float* spectrum, int kernel_width,int N)
{
  int k;
  float smoothing_tmp[N+1];
  float extended[N+2*kernel_width+1];
  float window[kernel_width*2 +1];
  float win_sum = 0.f;
  fft_window(window,kernel_width*2+1,0);//Hann window

  if (kernel_width == 0) return;

  //Copy data over the extended array to contemplate edge cases
  //Initialize smothingbins_tmp
  for (k = 0; k <= N; ++k)
  {
    extended[k+kernel_width] = spectrum[k];
    smoothing_tmp[k] = 0.f; //Initialize with zeros
  }

  for (k = 0; k <= kernel_width*2; ++k)
  {
    win_sum += window[k];
  }


  for (k = 0; k <= N; ++k)
  {
    for(int l = 0; l <= kernel_width*2; ++l)
    {
      smoothing_tmp[k] += window[l]*extended[k+l]/win_sum;
    }
  }

  for (k = 0; k <= N; ++k)
  {
    spectrum[k] = smoothing_tmp[k];
  }
}

//With quadratic coefficients
void
spectral_smoothing_SG_quad(float* spectrum, int kernel_width,int N)
{
  int k;
  float smoothing_tmp[N+1];
  memset(smoothing_tmp,0,sizeof(float)*(N+1));
  float extended[N+1+2*kernel_width];

  if (kernel_width < 2 || kernel_width > 12) return;

  //Copy data over the extended array to contemplate edge cases
  for (k = 0; k <= N; ++k)
  {
    extended[k+kernel_width] = spectrum[k];
  }

  for (k = 0; k <= N; ++k)
  {
    for(int l = 0; l <= kernel_width*2; ++l)
    {
      switch(kernel_width)
      {
        case 2:
          smoothing_tmp[k] += savgol_quad_5[l]*extended[l+k];
          break;
        case 3:
          smoothing_tmp[k] += savgol_quad_7[l]*extended[l+k];
          break;
        case 4:
          smoothing_tmp[k] += savgol_quad_9[l]*extended[l+k];
          break;
        case 5:
          smoothing_tmp[k] += savgol_quad_11[l]*extended[l+k];
          break;
        case 6:
          smoothing_tmp[k] += savgol_quad_13[l]*extended[l+k];
          break;
        case 7:
          smoothing_tmp[k] += savgol_quad_15[l]*extended[l+k];
          break;
        case 8:
          smoothing_tmp[k] += savgol_quad_17[l]*extended[l+k];
          break;
        case 9:
          smoothing_tmp[k] += savgol_quad_19[l]*extended[l+k];
          break;
        case 10:
          smoothing_tmp[k] += savgol_quad_21[l]*extended[l+k];
          break;
        case 11:
          smoothing_tmp[k] += savgol_quad_23[l]*extended[l+k];
          break;
        case 12:
          smoothing_tmp[k] += savgol_quad_25[l]*extended[l+k];
          break;
      }
    }
  }

  for (k = 0; k <= N; ++k)
  {
    spectrum[k] = smoothing_tmp[k];
  }
}

//With quadric coefficients
void
spectral_smoothing_SG_quart(float* spectrum, int kernel_width,int N)
{
  int k;
  float smoothing_tmp[N+1];
  float extended[N+2*kernel_width+1];

  if (kernel_width < 3 || kernel_width > 4) return;

  //Initialize smothingbins_tmp
  for (k = 0; k <= N; ++k)
  {
    extended[k+kernel_width] = spectrum[k];
    smoothing_tmp[k] = 0.f; //Initialize with zeros
  }

  for (k = 0; k <= N; ++k)
  {
    for(int l = 0; l < kernel_width*2; ++l)
    {
      switch(kernel_width)
      {
        case 3:
          smoothing_tmp[k] += savgol_quart_7[l]*extended[l+k];
          break;
        case 4:
          smoothing_tmp[k] += savgol_quart_9[l]*extended[l+k];
          break;
      }
    }
  }

  for (k = 0; k <= N; ++k)
  {
    spectrum[k] = smoothing_tmp[k];
  }
}
