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

#include "denoise_parameters.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define FROM_DB_TO_CV(gain_db) (expf(gain_db / 10.f * logf(10.f)))

struct DenoiseParameters
{
    float *enable;
    float *learn_noise;
    float *residual_listen;
    float *reduction_amount;
    float *release_time;
    float *masking_ceiling_limit;
    float *whitening_factor;
    float *transient_threshold;
    float *noise_rescale;
};

float get_plugin_parameters(DenoiseParameters *self, int parameter_type)
{
    switch ((ParameterType)parameter_type)
    {
    case REDUCTION_AMOUNT:
        return FROM_DB_TO_CV(-1.f * *self->reduction_amount);
        break;
    case NOISE_RESCALE:
        return *self->noise_rescale;
        break;
    case RELEASE:
        return *self->release_time;
        break;
    case MASKING:
        return *self->masking_ceiling_limit;
        break;
    case WHITENING_FACTOR:
        return *self->whitening_factor / 100.f;
        break;
    case LEARN_NOISE:
        return *self->learn_noise;
        break;
    case RESIDUAL_LISTEN:
        return *self->residual_listen;
        break;
    case TRANSIENT_PROTECTION:
        return *self->transient_threshold;
        break;
    case ENABLE:
        return *self->enable;
        break;
    }

    return 0;
}

void set_plugin_parameters(DenoiseParameters *self, float *value, int parameter_type)
{
    switch ((ParameterType)parameter_type)
    {
    case REDUCTION_AMOUNT:
        self->reduction_amount = value;
        break;
    case NOISE_RESCALE:
        self->noise_rescale = value;
        break;
    case RELEASE:
        self->release_time = value;
        break;
    case MASKING:
        self->masking_ceiling_limit = value;
        break;
    case WHITENING_FACTOR:
        self->whitening_factor = value;
        break;
    case LEARN_NOISE:
        self->learn_noise = value;
        break;
    case RESIDUAL_LISTEN:
        self->residual_listen = value;
        break;
    case TRANSIENT_PROTECTION:
        self->transient_threshold = value;
        break;
    case ENABLE:
        self->enable = value;
        break;
    }
}

void plugin_parameters_free(DenoiseParameters *self)
{
    free(self);
}

DenoiseParameters *plugin_parameters_initialize()
{
    DenoiseParameters *self = (DenoiseParameters *)calloc(1, sizeof(DenoiseParameters));

    return self;
}