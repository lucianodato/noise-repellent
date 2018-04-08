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
* \file spectrum_smoother.c
* \author Luciano Dato
* \brief Contains a spectrum smoother abstraction
*/

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

    //Ensemble related
    //Spectrum
    float *noise_spectrum;
    float *signal_spectrum;

    //smoothing related
    float *smoothed_spectrum;      //power spectrum to be smoothed
    float *smoothed_spectrum_prev; //previous frame smoothed power spectrum for envelopes
} FFTdenoiser;