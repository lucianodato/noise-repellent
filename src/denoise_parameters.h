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

#ifndef DENOISE_PARAMETERS_H
#define DENOISE_PARAMETERS_H

typedef enum
{
    REDUCTION_AMOUNT = 0,
    NOISE_RESCALE = 1,
    RELEASE = 2,
    MASKING = 3,
    WHITENING_FACTOR = 4,
    LEARN_NOISE = 5,
    RESIDUAL_LISTEN = 6,
    TRANSIENT_PROTECTION = 7,
    ENABLE = 8
} ParameterType;

typedef struct DenoiseParameters DenoiseParameters;

void plugin_parameters_free(DenoiseParameters *self);
float get_plugin_parameters(DenoiseParameters *self, int parameter_type);
void set_plugin_parameters(DenoiseParameters *self, float *value, int parameter_type);
DenoiseParameters *plugin_parameters_initialize();

#endif