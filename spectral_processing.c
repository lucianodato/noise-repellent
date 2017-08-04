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
#include "denoise_gain.c"
#include "masking.c"

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
  float fft_magnitude_db[fft_size_2+1];
  float peak_threshold_db = to_dB(SP_THRESH);
  int max_bin = MIN(freq_to_bin(SP_MAX_FREQ,samp_rate,fft_size_2),fft_size_2+1);
  int min_bin = MAX(freq_to_bin(SP_MIN_FREQ,samp_rate,fft_size_2),0);
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

//------------GAIN AND THRESHOLD CALCULATION---------------

void
preprocessing(float noise_thresholds_offset, float* noise_thresholds_scaled,
							float* smoothed_spectrum,	float* smoothed_spectrum_prev, int fft_size_2,	float* prev_beta, float* bark_z, float* absolute_thresholds, float* SSF,	float* max_masked, float* min_masked, float release_coeff)
{
	int k;

	//PREPROCESSING

	//------OVERSUBTRACTION------

	//CALCULATION OF ALPHA WITH MASKING THRESHOLDS USING VIRAGS METHOD
	// float alpha[fft_size_2+1];
	//
	// compute_alpha_and_beta(fft_p2, noise_thresholds_p2, fft_size_2, alpha, bark_z,
	// 											 absolute_thresholds, SSF, max_masked, min_masked);
	//Virag requires alphas to be smoothed over time and frequency? TODO


	//TODO
	// transient_preservation_coeff = transient_preservation(fft_p2, fft_p2_prev_tpres,
	// 																											fft_size_2);
	// memcpy(fft_p2_prev_tpres,fft_p2,sizeof(float)*(fft_size_2+1));


	//Scale noise thresholds (equals applying an oversubtraction factor in spectral subtraction)
	for (k = 0; k <= fft_size_2; k++)
	{
		noise_thresholds_scaled[k] *= noise_thresholds_offset;//* alpha[k] * transient_preservation_coeff;
	}

	//------SMOOTHING DETECTOR------

	/*Time smoothing between current and past Gk (similar effect to ephraim and malah)
		Here is done by applying a release envelope to signal power spectrum
		The best option here is to adaptively smooth 2D spectral components so it will require a biger buffer
		as suggested by Lukin in Suppression of Musical Noise Artifacts in Audio Noise Reduction by Adaptive 2D Filtering
	*/
	apply_envelope(smoothed_spectrum, smoothed_spectrum_prev, fft_size_2, release_coeff);

	// This adaptive method is based on SPECTRAL SUBTRACTION WITH ADAPTIVE AVERAGING OF THE GAIN FUNCTION
	// Not working correctly yet
	// spectrum_adaptive_time_smoothing(fft_size_2, smoothed_spectrum_prev, smoothed_spectrum,
	// 																 noise_thresholds_scaled, prev_beta, release_coeff);

	memcpy(smoothed_spectrum_prev,smoothed_spectrum,sizeof(float)*(fft_size_2+1));
}

void
spectral_gain(float* smoothed_spectrum, float* noise_thresholds_scaled, int fft_size_2,
							float adaptive, float* Gk)
{
	if(adaptive == 1.f)
	{
		power_subtraction(fft_size_2,	smoothed_spectrum, noise_thresholds_scaled, Gk);
	}
	else
	{
		spectral_gating(fft_size_2, smoothed_spectrum, noise_thresholds_scaled, Gk);
	}
}

void
postprocessing(int fft_size_2, int fft_size, float* fft_p2, float* output_fft_buffer,
							 float* input_fft_buffer_ps, float* input_fft_buffer_g,
							 float* output_fft_buffer_ps, float* output_fft_buffer_g,
							 fftwf_plan* forward_g, fftwf_plan* backward_g, fftwf_plan* forward_ps, 	float* Gk, float pf_threshold)
{

  int k;
  float postfilter[fft_size];

	//GAIN SMOOTHING USING A POSTFILTER
	//Compute the filter
	compute_post_filter(fft_size_2, fft_size, fft_p2, pf_threshold, postfilter, Gk);

	//Convolution using fft transform

	//Copy to fft buffers
	memcpy(input_fft_buffer_ps,postfilter,fft_size*sizeof(float));
	memcpy(input_fft_buffer_g,Gk,fft_size*sizeof(float));

	//FFT Analysis
	fftwf_execute(*forward_ps);
	fftwf_execute(*forward_g);

	//Multiply with the filter computed
	for (k = 0; k < fft_size; k++)
	{
    output_fft_buffer_g[k] *= output_fft_buffer_ps[k];
  }

	//FFT Synthesis (only gain needs to be synthetised)
	fftwf_execute(*backward_g);

	//Normalizing
	for (k = 0; k < fft_size; k++)
	{
		input_fft_buffer_g[k] = input_fft_buffer_g[k] / fft_size;
	}

	//Copy to orginal arrays
	memcpy(Gk,input_fft_buffer_g,fft_size*sizeof(float));

	///////////////////
}

void
denoised_calulation(int fft_size_2, int fft_size,	float* output_fft_buffer,
										float* denoised_spectrum, float* Gk)
{

  int k;

  //Apply the computed gain to the signal and store it in denoised array
  for (k = 0; k < fft_size; k++)
	{
    denoised_spectrum[k] = output_fft_buffer[k] * Gk[k];
  }
}

void
residual_calulation(int fft_size_2, int fft_size, float* output_fft_buffer,
										float* residual_spectrum, float* denoised_spectrum,
										float whitening_factor)
{

  int k;

  //Residual signal
  for (k = 0; k < fft_size; k++)
	{
   residual_spectrum[k] = output_fft_buffer[k] - denoised_spectrum[k];
  }

	////////////POSTPROCESSING RESIDUAL
	//Whitening (residual spectrum more similar to white noise)
	if(whitening_factor > 0.f)
	{
		whitening(residual_spectrum,whitening_factor,fft_size);
	}
	////////////
}

void
final_spectrum_ensemble(int fft_size_2, int fft_size, float* output_fft_buffer,
												float* residual_spectrum, float* denoised_spectrum,
												float reduction_amount, float wet_dry, float noise_listen)
{
  int k;

	//OUTPUT RESULTS using smooth bypass and parametric sustraction
	if (noise_listen == 0.f)
	{
	//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k < fft_size; k++)
		{
			output_fft_buffer[k] =  (1.f-wet_dry) * output_fft_buffer[k] + (denoised_spectrum[k] + residual_spectrum[k]*reduction_amount) * wet_dry;
		}
	}
	else
	{
		//Output noise only
		for (k = 0; k < fft_size; k++)
		{
			output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + residual_spectrum[k] * wet_dry;
		}
	}
}
