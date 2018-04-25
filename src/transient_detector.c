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
* \file transient_detector.c
* \author Luciano Dato
* \brief Contains a transient detector abstraction
*/

/**
* Gain estimation struct.
*/
typedef struct
{
    //General parameters
    int fft_size;
    int fft_size_2;

    float *spectrum;

    //Transient preservation related
    float *transient_preserv_prev; //previous frame smoothed power spectrum for envelopes
    float tp_r_mean;
    bool transient_present;
    float tp_window_count;
} Tdetector;

//------TRANSIENT PROTECTION------

/**
* Transient detection using a rolling mean thresholding over the spectral flux of
* the signal. Using more heuristics like high frequency content and others like the ones
* anylised by Dixon in 'Simple Spectrum-Based Onset Detection' would be better. Onset
* detection is explained thoroughly in 'A tutorial on onset detection in music signals' * by Bello.
*/
bool transient_detection(float *fft_p2, float *transient_preserv_prev, float fft_size_2,
                         float *tp_window_count, float *tp_r_mean, float transient_protection)
{
    float adapted_threshold, reduction_function;

    //Transient protection by forcing wiener filtering when an onset is detected
    reduction_function = spectral_flux(fft_p2, transient_preserv_prev, fft_size_2);
    //reduction_function = high_frequency_content(fft_p2, fft_size_2);

    *(tp_window_count) += 1.f;

    if (*(tp_window_count) > 1.f)
    {
        *(tp_r_mean) += ((reduction_function - *(tp_r_mean)) / *(tp_window_count));
    }
    else
    {
        *(tp_r_mean) = reduction_function;
    }

    adapted_threshold = (TP_UPPER_LIMIT - transient_protection) * *(tp_r_mean);

    memcpy(transient_preserv_prev, fft_p2, sizeof(float) * (fft_size_2 + 1));

    if (reduction_function > adapted_threshold)
    {
        return true;
    }
    else
    {
        return false;
    }
}