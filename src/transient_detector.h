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
* \file transient_detector.h
* \author Luciano Dato
* \brief Contains a transient detector abstraction
*/

#ifndef TRANSIENT_DETECTOR_H
#define TRANSIENT_DETECTOR_H

#include <stdbool.h>

/**
* Gain estimation struct.
*/
typedef struct TransientDetector TransientDetector;

float spectral_flux(float *spectrum, float *spectrum_prev, float N);
bool transient_detector_run(TransientDetector *self, float transient_threshold);
void transient_detector_reset(TransientDetector *self);
TransientDetector *transient_detector_initialize(int fft_size);
void transient_detector_free(TransientDetector *self);

#endif