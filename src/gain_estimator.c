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
* \file fft_denoiser.c
* \author Luciano Dato
* \brief Contains a the spectral noise reduction abstraction
*/

#include <float.h>
#include <math.h>

//General spectral subtraction configuration
#define GAMMA1 2.f
#define GAMMA2 0.5f

/**
* Gain estimation struct.
*/
typedef struct
{
    //General parameters
    int fft_size;
    int fft_size_2;
    int samp_rate;
    int hop;

    //Reduction gains
    float *Gk; //definitive gain

    //Ensemble related
    //Spectrum
    float *noise_spectrum;
    float *signal_spectrum;

	Mestimator masking_estimation;

    //smoothing related
    float *smoothed_spectrum;      //power spectrum to be smoothed
    float *smoothed_spectrum_prev; //previous frame smoothed power spectrum for envelopes

    //Transient preservation related
    float *transient_preserv_prev; //previous frame smoothed power spectrum for envelopes
    float tp_r_mean;
    bool transient_present;
    float tp_window_count;

    //whitening related
    float *residual_max_spectrum;
    float max_decay_rate;
    float whitening_window_count;

    //masking
    float *bark_z;
    float *absolute_thresholds; //absolute threshold of hearing
    float *SSF;
    float *spl_reference_values;
    float *unity_gain_bark_spectrum;
    float *spreaded_unity_gain_bark_spectrum;
    float *alpha_masking;
    float *beta_masking;
} FFTdenoiser;


/**
* Wiener substraction supression rule. Outputs the filter mirrored around nyquist.
* \param fft_size_2 is half of the fft size
* \param noise_thresholds is the threshold for each corresponding power spectum value
* /param spectrum is the power spectum array
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
*/
void wiener_subtraction(int fft_size_2, float *spectrum, float *noise_thresholds, float *Gk)
{
	int k;

	for (k = 0; k <= fft_size_2; k++)
	{
		if (noise_thresholds[k] > FLT_MIN)
		{
			if (spectrum[k] > noise_thresholds[k])
			{
				Gk[k] = (spectrum[k] - noise_thresholds[k]) / spectrum[k];
			}
			else
			{
				Gk[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		Gk[(2 * fft_size_2) - k] = Gk[k];
	}
}

/**
* Power substraction supression rule. Outputs the filter mirrored around nyquist.
* \param fft_size_2 is half of the fft size
* \param spectrum is the power spectum array
* \param noise_thresholds is the threshold for each corresponding power spectum value
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
*/
void power_subtraction(int fft_size_2, float *spectrum, float *noise_thresholds, float *Gk)
{
	int k;

	for (k = 0; k <= fft_size_2; k++)
	{
		if (noise_thresholds[k] > FLT_MIN)
		{
			if (spectrum[k] > noise_thresholds[k])
			{
				Gk[k] = sqrtf((spectrum[k] - noise_thresholds[k]) / spectrum[k]);
			}
			else
			{
				Gk[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		Gk[(2 * fft_size_2) - k] = Gk[k];
	}
}

/**
* Magnitude substraction supression rule. Outputs the filter mirrored around nyquist.
* \param fft_size_2 is half of the fft size
* \param spectrum is the power spectum array
* \param noise_thresholds is the threshold for each corresponding power spectum value
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
*/
void magnitude_subtraction(int fft_size_2, float *spectrum, float *noise_thresholds, float *Gk)
{
	int k;

	for (k = 0; k <= fft_size_2; k++)
	{
		if (noise_thresholds[k] > FLT_MIN)
		{
			if (spectrum[k] > noise_thresholds[k])
			{
				Gk[k] = (sqrtf(spectrum[k]) - sqrtf(noise_thresholds[k])) / sqrtf(spectrum[k]);
			}
			else
			{
				Gk[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		Gk[(2 * fft_size_2) - k] = Gk[k];
	}
}

/**
* Gating with hard knee supression rule. Outputs the filter mirrored around nyquist.
* \param fft_size_2 is half of the fft size
* \param spectrum is the power spectum array
* \param noise_thresholds is the threshold for each corresponding power spectum value
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
*/
void spectral_gating(int fft_size_2, float *spectrum, float *noise_thresholds, float *Gk)
{
	int k;

	for (k = 0; k <= fft_size_2; k++)
	{
		if (noise_thresholds[k] > FLT_MIN)
		{
			// //Without knee
			if (spectrum[k] >= noise_thresholds[k])
			{
				//over the threshold
				Gk[k] = 1.f;
			}
			else
			{
				//under the threshold
				Gk[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		Gk[(2 * fft_size_2) - k] = Gk[k];
	}
}

/**
* Generalized spectral subtraction supression rule. This version uses an array of alphas and betas. Outputs the filter mirrored around nyquist. GAMMA defines what type of spectral Subtraction is used. GAMMA1=GAMMA2=1 is magnitude substaction. GAMMA1=2 GAMMA2=0.5 is power Subtraction. GAMMA1=2 GAMMA2=1 is wiener filtering.
* \param fft_size_2 is half of the fft size
* \param alpha is the array of oversubtraction factors for each bin
* \param beta is the array of the spectral flooring factors for each bin
* \param spectrum is the power spectum array
* \param noise_thresholds is the threshold for each corresponding power spectum value
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
*/
void denoise_gain_gss(int fft_size_2, float *alpha, float *beta, float *spectrum,
					  float *noise_thresholds, float *Gk)
{
	int k;

	for (k = 0; k <= fft_size_2; k++)
	{
		if (spectrum[k] > FLT_MIN)
		{
			if (powf((noise_thresholds[k] / spectrum[k]), GAMMA1) < (1.f / (alpha[k] + beta[k])))
			{
				Gk[k] = MAX(powf(1.f - (alpha[k] * powf((noise_thresholds[k] / spectrum[k]), GAMMA1)), GAMMA2), 0.f);
			}
			else
			{
				Gk[k] = MAX(powf(beta[k] * powf((noise_thresholds[k] / spectrum[k]), GAMMA1), GAMMA2), 0.f);
			}
		}
		else
		{
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		Gk[(2 * fft_size_2) - k] = Gk[k];
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
void compute_alpha_and_beta(float *fft_p2, float *noise_thresholds_p2, int fft_size_2,
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