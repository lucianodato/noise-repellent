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
* \file gain_estimator.c
* \author Luciano Dato
* \brief Contains a the reduction gain estimator abstraction
*/

#include "masking_estimator.c"
#include "spectrum_smoother.c"
#include "transient_detector.c"
#include <float.h>
#include <math.h>
#include <stdbool.h>

//General spectral subtraction configuration
#define GAMMA1 2.f
#define GAMMA2 0.5f

//masking thresholds values recomended by virag
#define ALPHA_MAX 6.f
#define ALPHA_MIN 1.f
#define BETA_MAX 0.02f
#define BETA_MIN 0.0f

/**
* Gain estimation struct.
*/
typedef struct
{
	//General parameters
	int fft_size;
	int half_fft_size;
	int samp_rate;
	int hop;

	//Arrays
	float *gain_spectrum; //definitive gain
	float *noise_spectrum;
	float *signal_spectrum;
	float *alpha;
	float *beta;
	float *masking_thresholds;
	float *clean_signal_estimation;

	bool transient_detected;

	MaskingEstimator *masking_estimation;
	TransientDetector *transient_detection;
	// SpectralSmoother *spectrum_smoothing;
} GainEstimator;

/**
* Wiener substraction supression rule. Outputs the filter mirrored around nyquist.
*/
void wiener_subtraction(GainEstimator *self)
{
	int k;

	for (k = 0; k <= self->half_fft_size; k++)
	{
		if (self->noise_spectrum[k] > FLT_MIN)
		{
			if (self->signal_spectrum[k] > self->noise_spectrum[k])
			{
				self->gain_spectrum[k] = (self->signal_spectrum[k] - self->noise_spectrum[k]) / self->signal_spectrum[k];
			}
			else
			{
				self->gain_spectrum[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			self->gain_spectrum[k] = 1.f;
		}
	}
}

/**
* Power substraction supression rule. Outputs the filter mirrored around nyquist.
*/
void power_subtraction(GainEstimator *self)
{
	int k;

	for (k = 0; k <= self->half_fft_size; k++)
	{
		if (self->noise_spectrum[k] > FLT_MIN)
		{
			if (self->signal_spectrum[k] > self->noise_spectrum[k])
			{
				self->gain_spectrum[k] = sqrtf((self->signal_spectrum[k] - self->noise_spectrum[k]) / self->signal_spectrum[k]);
			}
			else
			{
				self->gain_spectrum[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			self->gain_spectrum[k] = 1.f;
		}
	}
}

/**
* Magnitude substraction supression rule. Outputs the filter mirrored around nyquist.
*/
void magnitude_subtraction(GainEstimator *self)
{
	int k;

	for (k = 0; k <= self->half_fft_size; k++)
	{
		if (self->noise_spectrum[k] > FLT_MIN)
		{
			if (self->signal_spectrum[k] > self->noise_spectrum[k])
			{
				self->gain_spectrum[k] = (sqrtf(self->signal_spectrum[k]) - sqrtf(self->noise_spectrum[k])) / sqrtf(self->signal_spectrum[k]);
			}
			else
			{
				self->gain_spectrum[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			self->gain_spectrum[k] = 1.f;
		}
	}
}

/**
* Gating with hard knee supression rule. Outputs the filter mirrored around nyquist.
*/
void spectral_gating(GainEstimator *self)
{
	int k;

	for (k = 0; k <= self->half_fft_size; k++)
	{
		if (self->noise_spectrum[k] > FLT_MIN)
		{
			// //Without knee
			if (self->signal_spectrum[k] >= self->noise_spectrum[k])
			{
				//over the threshold
				self->gain_spectrum[k] = 1.f;
			}
			else
			{
				//under the threshold
				self->gain_spectrum[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			self->gain_spectrum[k] = 1.f;
		}
	}
}

/**
* Generalized spectral subtraction supression rule. This version uses an array of alphas and betas. Outputs the filter mirrored around nyquist. 
* GAMMA defines what type of spectral Subtraction is used. 
* GAMMA1=GAMMA2=1 is magnitude substaction. 
* GAMMA1=2 GAMMA2=0.5 is power Subtraction. 
* GAMMA1=2 GAMMA2=1 is wiener filtering.
*/
void denoise_gain_gss(GainEstimator *self)
{
	int k;

	for (k = 0; k <= self->half_fft_size; k++)
	{
		if (self->signal_spectrum[k] > FLT_MIN)
		{
			if (powf((self->noise_spectrum[k] / self->signal_spectrum[k]), GAMMA1) < (1.f / (self->alpha[k] + self->beta[k])))
			{
				self->gain_spectrum[k] = MAX(powf(1.f - (self->alpha[k] * powf((self->noise_spectrum[k] / self->signal_spectrum[k]), GAMMA1)), GAMMA2), 0.f);
			}
			else
			{
				self->gain_spectrum[k] = MAX(powf(self->beta[k] * powf((self->noise_spectrum[k] / self->signal_spectrum[k]), GAMMA1), GAMMA2), 0.f);
			}
		}
		else
		{
			//Otherwise we keep everything as is
			self->gain_spectrum[k] = 1.f;
		}
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
*/
void compute_alpha_and_beta(GainEstimator *self, float masking_ceiling_limit, float masking_floor_limit)
{
	int k;
	float normalized_value;

	//Noise masking threshold must be computed from a clean signal
	//therefor we aproximate a clean signal using a power Sustraction over
	//the original noisy one

	//basic Power Sustraction to estimate clean signal
	for (k = 0; k <= self->half_fft_size; k++)
	{
		self->clean_signal_estimation[k] = MAX(self->signal_spectrum[k] - self->noise_spectrum[k], FLT_MIN);
	}

	//Now we can compute noise masking threshold
	//then we copy the masking thresholds values to this object masking threshold array
	compute_masking_thresholds(self->masking_estimation, self->signal_spectrum, self->masking_thresholds);

	//First we need the maximun and the minimun value of the masking threshold
	float max_masked_tmp = max_spectral_value(self->masking_thresholds, self->half_fft_size);
	float min_masked_tmp = min_spectral_value(self->masking_thresholds, self->half_fft_size);

	for (k = 0; k <= self->half_fft_size; k++)
	{
		//new alpha and beta vector
		if (self->masking_thresholds[k] == max_masked_tmp)
		{
			self->alpha[k] = ALPHA_MIN;
			self->beta[k] = BETA_MIN;
		}
		if (self->masking_thresholds[k] == min_masked_tmp)
		{
			self->alpha[k] = masking_ceiling_limit;
			self->beta[k] = masking_floor_limit;
		}
		if (self->masking_thresholds[k] < max_masked_tmp && self->masking_thresholds[k] > min_masked_tmp)
		{
			//Linear interpolation of the value between max and min masked threshold values
			normalized_value = (self->masking_thresholds[k] - min_masked_tmp) / (max_masked_tmp - min_masked_tmp);

			self->alpha[k] = (1.f - normalized_value) * ALPHA_MIN + normalized_value * masking_ceiling_limit;
			self->beta[k] = (1.f - normalized_value) * BETA_MIN + normalized_value * masking_floor_limit;
		}
	}
}

void g_e_run(GainEstimator *self, float *signal_spectrum, float *gain_spectrum, float transient_threshold,
			 float masking_ceiling_limit, float release, float noise_rescale)
{
	int k;

	memcpy(self->signal_spectrum, signal_spectrum, sizeof(float) * self->half_fft_size + 1);

	//PREPROCESSING - PRECALCULATIONS

	//------TRANSIENT DETECTION------

	if (transient_threshold > 1.f)
	{
		self->transient_detected = t_d_run(self->transient_detection, transient_threshold);
	}

	//COMPUTING OF ALPHA WITH MASKING THRESHOLDS USING VIRAGS METHOD
	if (masking_ceiling_limit > 1.f)
	{
		compute_alpha_and_beta(self, masking_ceiling_limit, 0.f); //temporary value to imitate current state
	}
	else
	{
		initialize_array(self->alpha, 1.f, self->half_fft_size + 1); //This avoids incorrect results when moving sliders rapidly
	}

	//------OVERSUBTRACTION------

	//Scale noise thresholds (equals applying an oversubtraction factor in spectral subtraction)
	for (k = 0; k <= self->half_fft_size; k++)
	{
		self->noise_spectrum[k] = self->noise_spectrum[k] * noise_rescale * self->alpha[k];
	}

	//------SMOOTHING DETECTOR------

	// s_s_run(self->spectrum_smoothing, release);

	//------REDUCTION GAINS------

	//Protect transient by avoiding smoothing if present
	if (self->transient_detected && transient_threshold > 1.f)
	{
		wiener_subtraction(self);
	}
	else
	{
		spectral_gating(self);
	}

	//Copy the obtained values and mirror the gain array
	memcpy(gain_spectrum, self->gain_spectrum, sizeof(float) * self->half_fft_size + 1);
	for (k = 1; k < self->half_fft_size; k++)
	{
		gain_spectrum[(2 * self->half_fft_size) - k] = self->gain_spectrum[k];
	}
}

/**
* Reset dynamic arrays to zero.
*/
void g_e_reset(GainEstimator *self)
{
	//Reset all arrays
	initialize_array(self->signal_spectrum, 0.f, self->half_fft_size + 1);
	initialize_array(self->noise_spectrum, 0.f, self->half_fft_size + 1);
	initialize_array(self->gain_spectrum, 1.f, self->half_fft_size + 1);
	initialize_array(self->alpha, 1.f, self->half_fft_size + 1);
	initialize_array(self->beta, 0.f, self->half_fft_size + 1);
	initialize_array(self->masking_thresholds, 0.f, self->half_fft_size + 1);
	initialize_array(self->clean_signal_estimation, 0.f, self->half_fft_size + 1);
}

/**
* Gain estimator initialization and configuration.
*/
GainEstimator *
g_e_init(int fft_size, int samp_rate, int hop)
{
	//Allocate object
	GainEstimator *self = (GainEstimator *)malloc(sizeof(GainEstimator));

	//Configuration
	self->fft_size = fft_size;
	self->half_fft_size = self->fft_size / 2;
	self->samp_rate = samp_rate;
	self->hop = hop;

	//spectrum allocation
	self->signal_spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->gain_spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->noise_spectrum = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->alpha = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->beta = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->masking_thresholds = (float *)calloc((self->half_fft_size + 1), sizeof(float));
	self->clean_signal_estimation = (float *)calloc((self->half_fft_size + 1), sizeof(float));

	//Reset all values
	g_e_reset(self);

	self->masking_estimation = m_e_init(self->fft_size, self->samp_rate);
	self->transient_detection = t_d_init(self->fft_size);
	// self->spectrum_smoothing = s_s_init();

	return self;
}

/**
* Free allocated memory.
*/
void g_e_free(GainEstimator *self)
{
	free(self->noise_spectrum);
	free(self->gain_spectrum);
	free(self->signal_spectrum);
	free(self->alpha);
	free(self->beta);
	free(self->masking_thresholds);
	free(self->clean_signal_estimation);
	m_e_free(self->masking_estimation);
	t_d_free(self->transient_detection);
	// s_s_free(self->spectrum_smoothing);
	free(self);
}