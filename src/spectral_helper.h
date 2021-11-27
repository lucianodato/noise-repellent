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
* \file spectral_helper.c
* \author Luciano Dato
* \brief Extra methods used by others. This keeps clean other files.
*/

#ifndef SPECTRAL_HELPER_H
#define SPECTRAL_HELPER_H

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

/**
* To initialize an array to a single value in all positions.
* \param array the array to initialize
* \param value the value to copy to every position in the array
* \param size the size of the array
*/
static void initialize_spectrum(float *array, float value, int size)
{
	for (int k = 0; k < size; k++)
	{
		array[k] = value;
	}
}

/**
* blackman window values computing.
* \param k bin number
* \param N fft size
*/
static float blackman(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.42 - 0.5 * cosf(2.f * M_PI * p) + 0.08 * cosf(4.f * M_PI * p);
}

/**
* hanning window values computing.
* \param k bin number
* \param N fft size
*/
static float hanning(int k, int N)
{
	float p = ((float)(k)) / ((float)(N));
	return 0.5 - 0.5 * cosf(2.f * M_PI * p);
}

/**
* hamming window values computing.
* \param k bin number
* \param N fft size
*/
static float hamming(int k, int N)
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
static float vorbis(int k, int N)
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
static void fft_window(float *window, int N, int window_type)
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

#endif
