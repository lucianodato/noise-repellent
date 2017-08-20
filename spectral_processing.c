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

#define ONSET_THRESH 80.f //For onset detection

//------------GAIN AND THRESHOLD CALCULATION---------------

void
preprocessing(float noise_thresholds_offset, float* fft_p2,
							float* noise_thresholds_p2,
							float* noise_thresholds_scaled, float* smoothed_spectrum,
							float* smoothed_spectrum_prev, int fft_size_2,
							float* prev_beta, float* bark_z, float* absolute_thresholds, float* SSF,
							float release_coeff, float* spreaded_unity_gain_bark_spectrum,
							float* spl_reference_values, float* alpha_masking, float* beta_masking,
							float masking_value, float adaptive_state, float reduction_value)
{
	int k;

	//PREPROCESSING

	//CALCULATION OF ALPHA WITH MASKING THRESHOLDS USING VIRAGS METHOD

	if(masking_value > 1.f && adaptive_state == 0.f){ //Only when adaptive is off
		compute_alpha_and_beta(fft_p2, noise_thresholds_p2, fft_size_2,
													 alpha_masking, beta_masking, bark_z, absolute_thresholds,
													 SSF, spreaded_unity_gain_bark_spectrum, spl_reference_values,
												 	 masking_value, reduction_value);
	}

	//------OVERSUBTRACTION------

	//Scale noise thresholds (equals applying an oversubtraction factor in spectral subtraction)
	for (k = 0; k <= fft_size_2; k++)
	{
		if(adaptive_state == 0.f)
		{
			noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset * alpha_masking[k];
		}
		else
		{
			noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset;
		}
	}

	//------SMOOTHING DETECTOR------

	if(adaptive_state == 0.f) //Only when adaptive is off
	{
		memcpy(smoothed_spectrum,fft_p2,sizeof(float)*(fft_size_2+1));

		apply_time_envelope(smoothed_spectrum, smoothed_spectrum_prev, fft_size_2, release_coeff);

		// This adaptive method is based on SPECTRAL SUBTRACTION WITH ADAPTIVE AVERAGING OF THE GAIN FUNCTION
		// spectrum_adaptive_time_smoothing(fft_size_2, smoothed_spectrum_prev, smoothed_spectrum,
		// 																 noise_thresholds_scaled, prev_beta, 1.f-release_coeff);

		memcpy(smoothed_spectrum_prev,smoothed_spectrum,sizeof(float)*(fft_size_2+1));
	}
}

void
spectral_gain(float* fft_p2, float* noise_thresholds_p2, float* noise_thresholds_scaled,
							float* smoothed_spectrum, int fft_size_2, float adaptive, float* Gk,
							float* transient_preserv_prev, float* reduction_function_prev,
							float* tp_window_count, float* tp_r_mean, float transient_protection)
{
	//------TRANSIENT PROTECTION------
	float adapted_threshold, reduction_function;

	//Transient protection by forcing wiener filtering when an onset is detected
	reduction_function = spectral_flux(fft_p2, transient_preserv_prev, fft_size_2);
	//reduction_function = high_frequency_content(fft_p2, fft_size_2);

	//adaptive thresholding (using rolling mean)
	*(tp_window_count) += 1.f;

	if(*(tp_window_count) > 1.f)
	{
		*(tp_r_mean) += ((reduction_function - *(reduction_function_prev))/ *(tp_window_count));
	}
	else
	{
		*(tp_r_mean) = reduction_function;
	}

	adapted_threshold = ONSET_THRESH + *(tp_r_mean);

	*(reduction_function_prev) = *(tp_r_mean);
	memcpy(transient_preserv_prev,fft_p2,sizeof(float)*(fft_size_2+1));

	//------REDUCTION GAINS------

	//Get reduction to apply
	if(adaptive == 1.f)
	{
		power_subtraction(fft_size_2, fft_p2, noise_thresholds_scaled, Gk);
	}
	else
	{
		//Protect transient by avoiding smoothing if present
		if (reduction_function > adapted_threshold && transient_protection == 1.f)
		{
			spectral_gating(fft_size_2, fft_p2, noise_thresholds_scaled, Gk);

			printf("%f", reduction_function);
			printf("%s", "   ");
			printf("%f", adapted_threshold);
			printf("%s", "   ");
			printf("%f\n", *(tp_r_mean));
		}
		else
		{
			spectral_gating(fft_size_2, smoothed_spectrum, noise_thresholds_scaled, Gk);
		}
	}
}

void
denoised_calulation(int fft_size,	float* output_fft_buffer,
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
residual_calulation(int fft_size, float* output_fft_buffer,
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
		spectral_whitening(residual_spectrum,whitening_factor,fft_size);
	}
	////////////
}

void
final_spectrum_ensemble(int fft_size, float* final_spectrum,
												float* residual_spectrum, float* denoised_spectrum,
												float reduction_amount, float noise_listen)
{
  int k;

	//OUTPUT RESULTS using smooth bypass and parametric sustraction
	if (noise_listen == 0.f)
	{
	//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k < fft_size; k++)
		{
			final_spectrum[k] = denoised_spectrum[k] + residual_spectrum[k]*reduction_amount;
		}
	}
	else
	{
		//Output noise only
		for (k = 0; k < fft_size; k++)
		{
			final_spectrum[k] = residual_spectrum[k];
		}
	}
}

void
soft_bypass(float* final_spectrum, float* output_fft_buffer, float wet_dry,
                 int fft_size)
{
  int k;

  for (k = 0; k < fft_size; k++)
  {
    output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + final_spectrum[k] * wet_dry;
  }
}
