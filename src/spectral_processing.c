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
* \file spectral_processing.c
* \author Luciano Dato
* \brief All methods related to spectral processing the spectrum
*/

#include <float.h>
#include <math.h>

#include "estimate_noise_spectrum.c"
#include "denoise_gain.c"
#include "masking.c"

//------------GAIN AND THRESHOLD CALCULATION---------------

/**
* Includes every preprocessing or precomputing before the supression rule.
* \param noise_thresholds_offset the scaling of the thresholds setted by the user
* \param fft_p2 the power spectrum of current frame
* \param noise_thresholds_p2 the noise thresholds for each bin estimated previously
* \param noise_thresholds_scaled the noise thresholds for each bin estimated previously scaled by the user
* \param smoothed_spectrum current power specturm with time smoothing applied
* \param smoothed_spectrum_prev the power specturm with time smoothing applied of previous frame
* \param fft_size_2 is half of the fft size
* \param prev_beta beta of previous frame for adaptive smoothing (not used yet)
* \param bark_z defines the bark to linear mapping for current spectrum config
* \param absolute_thresholds defines the absolute thresholds of hearing for current spectrum config
* \param SSF defines the spreading function matrix
* \param release_coeff release coefficient for time smoothing
* \param spreaded_unity_gain_bark_spectrum correction to be applied to SSF convolution
* \param spl_reference_values defines the reference values for each bin to convert from db to db SPL
* \param alpha_masking is the array of oversubtraction factors for each bin
* \param beta_masking is the array of the spectral flooring factors for each bin
* \param masking_value is the limit max oversubtraction to be computed
* \param adaptive flag that indicates if the noise is being estimated adaptively
* \param reduction_value is the limit max the spectral flooring to be computed
* \param transient_preserv_prev is the previous frame for spectral flux computing
* \param tp_window_count is the frame counter for the rolling mean thresholding for onset detection
* \param tp_r_mean is the rolling mean value for onset detection
* \param transient_present indicates if current frame is an onset or not (contains a transient)
* \param transient_protection is the flag that indicates whether transient protection is active or not
*/
void preprocessing(float noise_thresholds_offset, float *fft_p2,
				   float *noise_thresholds_p2,
				   float *noise_thresholds_scaled, float *smoothed_spectrum,
				   float *smoothed_spectrum_prev, int fft_size_2,
				   float *bark_z, float *absolute_thresholds, float *SSF,
				   float release_coeff, float *spreaded_unity_gain_bark_spectrum,
				   float *spl_reference_values, float *alpha_masking, float *beta_masking,
				   float masking_value, float adaptive_state, float reduction_value,
				   float *transient_preserv_prev, float *tp_window_count, float *tp_r_mean,
				   bool *transient_present, float transient_protection)
{
	int k;

	//PREPROCESSING - PRECALCULATIONS

	//------TRANSIENT DETECTION------

	if (transient_protection > 1.f)
	{
		*(transient_present) = transient_detection(fft_p2, transient_preserv_prev,
												   fft_size_2, tp_window_count, tp_r_mean,
												   transient_protection);
	}

	//CALCULATION OF ALPHA WITH MASKING THRESHOLDS USING VIRAGS METHOD

	if (masking_value > 1.f && adaptive_state == 0.f)
	{ //Only when adaptive is off
		compute_alpha_and_beta(fft_p2, noise_thresholds_p2, fft_size_2,
							   alpha_masking, beta_masking, bark_z, absolute_thresholds,
							   SSF, spreaded_unity_gain_bark_spectrum, spl_reference_values,
							   masking_value, reduction_value);
	}
	else
	{
		initialize_array(alpha_masking, 1.f, fft_size_2 + 1); //This avoids incorrect results when moving sliders rapidly
	}

	//------OVERSUBTRACTION------

	//Scale noise thresholds (equals applying an oversubtraction factor in spectral subtraction)
	for (k = 0; k <= fft_size_2; k++)
	{
		if (adaptive_state == 0.f)
		{
			noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset * alpha_masking[k];
		}
		else
		{
			noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset;
		}
	}

	//------SMOOTHING DETECTOR------

	if (adaptive_state == 0.f) //Only when adaptive is off
	{
		memcpy(smoothed_spectrum, fft_p2, sizeof(float) * (fft_size_2 + 1));

		apply_time_envelope(smoothed_spectrum, smoothed_spectrum_prev, fft_size_2, release_coeff);

		memcpy(smoothed_spectrum_prev, smoothed_spectrum, sizeof(float) * (fft_size_2 + 1));
	}
}

/**
* Computes the supression filter based on pre-processing data.
* \param fft_p2 the power spectrum of current frame
* \param noise_thresholds_p2 the noise thresholds for each bin estimated previously
* \param noise_thresholds_scaled the noise thresholds for each bin estimated previously scaled by the user
* \param smoothed_spectrum current power specturm with time smoothing applied
* \param fft_size_2 is half of the fft size
* \param adaptive flag that indicates if the noise is being estimated adaptively
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
* \param transient_protection is the flag that indicates whether transient protection is active or not
* \param transient_present indicates if current frame is an onset or not (contains a transient)
*/
void spectral_gain(float *fft_p2, float *noise_thresholds_p2, float *noise_thresholds_scaled,
				   float *smoothed_spectrum, int fft_size_2, float adaptive, float *Gk,
				   float transient_protection, bool transient_present)
{
	//------REDUCTION GAINS------

	//Get reduction to apply
	if (adaptive == 1.f)
	{
		power_subtraction(fft_size_2, fft_p2, noise_thresholds_scaled, Gk);
	}
	else
	{
		//Protect transient by avoiding smoothing if present
		if (transient_present && transient_protection > 1.f)
		{
			wiener_subtraction(fft_size_2, fft_p2, noise_thresholds_scaled, Gk);
		}
		else
		{
			spectral_gating(fft_size_2, smoothed_spectrum, noise_thresholds_scaled, Gk);
		}
	}
}

/**
* Applies the filter to the complex spectrum and gets the clean signal.
* \param fft_size size of the fft
* \param output_fft_buffer the unprocessed spectrum remaining in the fft buffer
* \param denoised_spectrum the spectrum of the cleaned signal
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
*/
void denoised_calulation(int fft_size, float *output_fft_buffer,
						 float *denoised_spectrum, float *Gk)
{
	int k;

	//Apply the computed gain to the signal and store it in denoised array
	for (k = 0; k < fft_size; k++)
	{
		denoised_spectrum[k] = output_fft_buffer[k] * Gk[k];
	}
}

/**
* Gets the residual signal of the reduction.
* \param fft_size size of the fft
* \param output_fft_buffer the unprocessed spectrum remaining in the fft buffer
* \param denoised_spectrum the spectrum of the cleaned signal
* \param whitening_factor the mix coefficient between whitened and not whitened residual spectrum
* \param residual_max_spectrum contains the maximun temporal value in each residual bin
* \param whitening_window_count counts frames to distinguish the first from the others
* \param max_decay_rate coefficient that sets the memory for each temporal maximun
*/
void residual_calulation(int fft_size, float *output_fft_buffer,
						 float *residual_spectrum, float *denoised_spectrum,
						 float whitening_factor, float *residual_max_spectrum,
						 float *whitening_window_count, float max_decay_rate)
{

	int k;

	//Residual signal
	for (k = 0; k < fft_size; k++)
	{
		residual_spectrum[k] = output_fft_buffer[k] - denoised_spectrum[k];
	}

	////////////POSTPROCESSING RESIDUAL
	//Whitening (residual spectrum more similar to white noise)
	if (whitening_factor > 0.f)
	{
		spectral_whitening(residual_spectrum, whitening_factor, fft_size,
						   residual_max_spectrum, whitening_window_count, max_decay_rate);
	}
	////////////
}

/**
* Mixes the cleaned signal with the residual taking into account the reduction configured
* by the user. Outputs the final signal or the residual only.
* \param fft_size size of the fft
* \param final_spectrum the spectrum to output from the plugin
* \param residual_spectrum the spectrum of the reduction residual
* \param denoised_spectrum the spectrum of the cleaned signal
* \param reduction_amount the amount of dB power to reduce setted by the user
* \param noise_listen control variable that decides whether to output the mixed noise reduced signal or the residual only
*/
void final_spectrum_ensemble(int fft_size, float *final_spectrum,
							 float *residual_spectrum, float *denoised_spectrum,
							 float reduction_amount, float noise_listen)
{
	int k;

	//OUTPUT RESULTS using smooth bypass and parametric subtraction
	if (noise_listen == 0.f)
	{
		//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k < fft_size; k++)
		{
			final_spectrum[k] = denoised_spectrum[k] + residual_spectrum[k] * reduction_amount;
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

/**
* Mixes unprocessed and processed signal to bypass softly.
* \param final_spectrum the spectrum to output from the plugin
* \param output_fft_buffer the unprocessed spectrum remaining in the fft buffer
* \param wet_dry mixing coefficient
* \param fft_size size of the fft
*/
void soft_bypass(float *final_spectrum, float *output_fft_buffer, float wet_dry,
				 int fft_size)
{
	int k;

	for (k = 0; k < fft_size; k++)
	{
		output_fft_buffer[k] = (1.f - wet_dry) * output_fft_buffer[k] + final_spectrum[k] * wet_dry;
	}
}
