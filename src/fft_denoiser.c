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
* \brief Contains an abstraction for a single fft spectrum denoising
*/

#ifndef FFT_DENOISER_C
#define FFT_DENOISER_C

#define WHITENING_DECAY_RATE 1000.f //Deacay in ms for max spectrum for whitening
#define WHITENING_FLOOR 0.02f		//Minumum max value posible

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#include "gain_estimator.c"
#include "noise_estimator.c"

/**
* FFT processor struct.
*/
typedef struct
{
	//General parameters
	int fft_size;
	int half_fft_size;
	int samp_rate;
	int hop;

	//Ensemble related
	float *fft_spectrum;
	float *processed_fft_spectrum;

	//Soft bypass
	float tau;			  //time constant for soft bypass
	float wet_dry_target; //softbypass target for softbypass
	float wet_dry;		  //softbypass coeff

	//Spectrum information arrays
	float *power_spectrum;
	float *phase_spectrum;
	float *magnitude_spectrum;

	float *gain_spectrum; //definitive reduction gain
	float *residual_spectrum;
	float *denoised_spectrum;
	float *whitened_residual_spectrum;

	GainEstimator *gain_estimation;
	NoiseEstimator *noise_estimation;

	//whitening related
	float *residual_max_spectrum;
	float max_decay_rate;
	float whitening_window_count;
} FFTDenoiser;

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
static void get_info_from_bins(float *fft_p2, float *fft_magnitude, float *fft_phase,
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
		fft_phase[k] = phase;	//This is not used but part of the STFT transform for generic use
	}
}

/**
* Verifies if the spectrum is full of zeros.
* \param spectrum the array to check
* \param N the size of the array (half the fft size plus 1)
*/
static bool is_empty(float *spectrum, int N)
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
* Updates the wet/dry mixing coefficient.
*/
void fft_processor_update_wetdry_target(FFTDenoiser *self, bool enable)
{
	//Softbypass targets in case of disabled or enabled
	if (enable)
	{ //if enabled
		self->wet_dry_target = 1.f;
	}
	else
	{ //if disabled
		self->wet_dry_target = 0.f;
	}
	//Interpolate parameters over time softly to bypass without clicks or pops
	self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry) + FLT_MIN;
}

/**
* Mixes unprocessed and processed signal to bypass softly.
*/
void fft_processor_soft_bypass(FFTDenoiser *self)
{
	int k;

	for (k = 0; k < self->fft_size; k++)
	{
		self->processed_fft_spectrum[k] = (1.f - self->wet_dry) * self->fft_spectrum[k] + self->processed_fft_spectrum[k] * self->wet_dry;
	}
}

/**
* Whitens the spectrum adaptively as proposed in 'Adaptive whitening for improved
* real-time audio onset detection' by Stowell and Plumbley. The idea here is that when
* residual noise resembles white noise the ear is able to precieve it as not so annoying.
* It uses a temporal max value for each bin and a decay factor as the memory regulator of
* that maximun value.
*/
void residual_spectrum_whitening(FFTDenoiser *self, float whitening_factor)
{
	self->whitening_window_count++;

	for (int k = 0; k < self->fft_size; k++)
	{
		if (self->whitening_window_count > 1.f)
		{
			self->residual_max_spectrum[k] = MAX(MAX(self->residual_spectrum[k], WHITENING_FLOOR), self->residual_max_spectrum[k] * self->max_decay_rate);
		}
		else
		{
			self->residual_max_spectrum[k] = MAX(self->residual_spectrum[k], WHITENING_FLOOR);
		}
	}

	for (int k = 0; k < self->fft_size; k++)
	{
		if (self->residual_spectrum[k] > FLT_MIN)
		{
			//Get whitened spectrum
			self->whitened_residual_spectrum[k] = self->residual_spectrum[k] / self->residual_max_spectrum[k];

			//Interpolate between whitened and non whitened residual
			self->residual_spectrum[k] = (1.f - whitening_factor) * self->residual_spectrum[k] + whitening_factor * self->whitened_residual_spectrum[k];
		}
	}
}

/**
* Applies the filter to the complex spectrum and gets the clean signal.
*/
void get_denoised_spectrum(FFTDenoiser *self)
{
	int k;

	//Apply the computed gain to the signal and store it in denoised array
	for (k = 0; k < self->fft_size; k++)
	{
		self->denoised_spectrum[k] = self->fft_spectrum[k] * self->gain_spectrum[k];
	}
}

/**
* Gets the residual signal of the reduction.
*/
void get_residual_spectrum(FFTDenoiser *self, float whitening_factor)
{
	int k;

	//Residual signal
	for (k = 0; k < self->fft_size; k++)
	{
		self->residual_spectrum[k] = self->fft_spectrum[k] - self->denoised_spectrum[k];
	}

	//Whitening (residual spectrum more similar to white noise)
	if (whitening_factor > 0.f)
	{
		residual_spectrum_whitening(self, whitening_factor);
	}
}

/**
* Mixes the cleaned signal with the residual taking into account the reduction configured
* by the user. Outputs the final signal or the residual only.
*/
void get_final_spectrum(FFTDenoiser *self, bool residual_listen, float reduction_amount)
{
	int k;

	//OUTPUT RESULTS using smooth bypass and parametric subtraction
	if (residual_listen)
	{
		//Output noise only
		for (k = 0; k < self->fft_size; k++)
		{
			self->processed_fft_spectrum[k] = self->residual_spectrum[k];
		}
	}
	else
	{
		//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k < self->fft_size; k++)
		{
			self->processed_fft_spectrum[k] = self->denoised_spectrum[k] +
											  self->residual_spectrum[k] * reduction_amount;
		}
	}
}

/**
* Runs the fft processing for current block.
*/
void fft_processor_run(FFTDenoiser *self, float *fft_spectrum, int enable, bool learn_noise, float whitening_factor,
					   float reduction_amount, bool residual_listen, float transient_threshold,
					   float masking_ceiling_limit, float release, float noise_rescale)
{
	fft_processor_update_wetdry_target(self, enable);

	memcpy(self->fft_spectrum, fft_spectrum, sizeof(float) * self->fft_size);

	//First get the power, magnitude and phase spectrum
	get_info_from_bins(self->power_spectrum, self->magnitude_spectrum,
					   self->phase_spectrum, self->half_fft_size,
					   self->fft_size, self->fft_spectrum);

	//If the spectrum is not silence
	if (!is_empty(self->power_spectrum, self->half_fft_size))
	{
		/*If selected estimate noise spectrum is based on selected portion of signal
        *do not process the signal
        */
		if (learn_noise)
		{
			//LEARN NOISE (Using power spectrum)
			noise_estimation_run(self->noise_estimation, self->power_spectrum);
		}
		else
		{
			//REDUCE NOISE OR LISTEN TO THE RESIDUAL
			if (is_noise_estimation_available(self->noise_estimation))
			{
				gain_estimation_run(self->gain_estimation, self->power_spectrum, self->gain_spectrum, transient_threshold,
									masking_ceiling_limit, release, noise_rescale);

				get_denoised_spectrum(self);

				get_residual_spectrum(self, whitening_factor);

				get_final_spectrum(self, residual_listen, reduction_amount);
			}
		}
	}

	//If bypassed mix unprocessed and processed signal softly
	fft_processor_soft_bypass(self);

	//Copy the processed spectrum to fft_spectrum
	memcpy(fft_spectrum, self->processed_fft_spectrum, sizeof(float) * self->fft_size);
}

/**
* Reset dynamic arrays to zero.
*/
void fft_processor_reset(FFTDenoiser *self)
{
	//Reset all arrays
	memset(self->fft_spectrum, 0.f, self->fft_size);
	memset(self->processed_fft_spectrum, 0.f, self->fft_size);
	memset(self->gain_spectrum, 1.f, self->fft_size);

	memset(self->power_spectrum, 0.f, self->half_fft_size + 1);
	memset(self->magnitude_spectrum, 0.f, self->half_fft_size + 1);
	memset(self->phase_spectrum, 0.f, self->half_fft_size + 1);

	memset(self->residual_max_spectrum, 0.f, self->fft_size);
	memset(self->denoised_spectrum, 0.f, self->fft_size);
	memset(self->residual_spectrum, 0.f, self->fft_size);
	memset(self->whitened_residual_spectrum, 0.f, self->fft_size);
	memset(self->gain_spectrum, 0.f, self->fft_size);

	self->whitening_window_count = 0.f;
}

/**
* FFT processor initialization and configuration.
*/
FFTDenoiser *
fft_processor_initialize(int fft_size, int samp_rate, int hop)
{
	//Allocate object
	FFTDenoiser *self = (FFTDenoiser *)malloc(sizeof(FFTDenoiser));

	//Configuration
	self->fft_size = fft_size;
	self->half_fft_size = self->fft_size / 2;
	self->samp_rate = samp_rate;
	self->hop = hop;

	//spectrum allocation
	self->fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->processed_fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));

	//soft bypass
	self->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->samp_rate));
	self->wet_dry = 0.f;

	//whitening related
	self->residual_max_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->max_decay_rate = expf(-1000.f / (((WHITENING_DECAY_RATE)*self->samp_rate) / self->hop));

	//final ensemble related
	self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->gain_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->whitened_residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));

	//Arrays for getting bins info
	self->power_spectrum = (float *)malloc((self->half_fft_size + 1) * sizeof(float));
	self->magnitude_spectrum = (float *)malloc((self->half_fft_size + 1) * sizeof(float));
	self->phase_spectrum = (float *)malloc((self->half_fft_size + 1) * sizeof(float));

	//Reset all values
	fft_processor_reset(self);

	//Noise estimator related
	self->noise_estimation = noise_estimation_initialize(self->fft_size);

	return self;
}

/**
* Free allocated memory.
*/
void fft_processor_free(FFTDenoiser *self)
{
	free(self->fft_spectrum);
	free(self->processed_fft_spectrum);
	free(self->power_spectrum);
	free(self->magnitude_spectrum);
	free(self->phase_spectrum);
	free(self->gain_spectrum);
	free(self->residual_spectrum);
	free(self->whitened_residual_spectrum);
	free(self->denoised_spectrum);
	free(self->residual_max_spectrum);
	noise_estimation_free(self->noise_estimation);
	free(self);
}

#endif