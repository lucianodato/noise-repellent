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
