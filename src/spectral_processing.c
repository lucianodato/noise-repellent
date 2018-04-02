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
* \file spectral_processing.c
* \author Luciano Dato
* \brief Contains an spectral processing abstraction
* \ in this case noise reduction
*/

#include "extra_functions.c"

/**
* Struct for the processing unit.
*/
typedef struct
{
    //General parameters
    int fft_size;
    int fft_size_2;
    int samp_rate;
    int hop;

    //Ensemble related
    float *original_fft_spectrum;
    float *processed_fft_spectrum;

    //Spectrum related
    float *power_spectrum;
    float *phase_spectrum;
    float *magnitude_spectrum;

    //Soft bypass
    float tau;            //time constant for soft bypass
    float wet_dry_target; //softbypass target for softbypass
    float wet_dry;        //softbypass coeff

    //Sdenoiser *denoiser;
} Sprocessor;

void sp_configure(Sprocessor *self, int fft_size, int samp_rate, int hop)
{
    self->fft_size = fft_size;
    self->fft_size_2 = self->fft_size / 2;
    self->samp_rate = samp_rate;
    self->hop = hop;
}

void sp_reset(Sprocessor *self)
{
    //Reset all arrays
    initialize_array(self->original_fft_spectrum, 0.f, self->fft_size);
    initialize_array(self->processed_fft_spectrum, 0.f, self->fft_size);
    initialize_array(self->power_spectrum, 0.f, self->fft_size_2 + 1);
    initialize_array(self->magnitude_spectrum, 0.f, self->fft_size_2 + 1);
    initialize_array(self->phase_spectrum, 0.f, self->fft_size_2 + 1);
}

Sprocessor *
sp_init(int fft_size, int samp_rate, int hop)
{
    //Allocate object
    Sprocessor *self = (Sprocessor *)malloc(sizeof(Sprocessor));

    sp_configure(self, fft_size, samp_rate, hop);

    //spectrum allocation
    self->original_fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));
    self->processed_fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));

    //Arrays for getting bins info
    self->power_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));
    self->magnitude_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));
    self->phase_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));

    //soft bypass
    self->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->samp_rate));
    self->wet_dry = 0.f;

    //Initialize arrays with zeros
    sp_reset(self);

    //Estimation initialization
    //sd_init(self->denoiser);

    return self;
}

void sp_free(Sprocessor *self)
{
    free(self->power_spectrum);
    free(self->magnitude_spectrum);
    free(self->phase_spectrum);
    free(self->original_fft_spectrum);
    free(self->processed_fft_spectrum);
    //sd_free(self->estimation);
    free(self);
}

void update_wetdry_target(Sprocessor *self, float* enable)
{
    //Softbypass targets in case of disabled or enabled
    if (*(enable) == 0.f)
    { //if disabled
        self->wet_dry_target = 0.f;
    }
    else
    { //if enabled
        self->wet_dry_target = 1.f;
    }
    //Interpolate parameters over time softly to bypass without clicks or pops
    self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry) + FLT_MIN;
}

/**
* Mixes unprocessed and processed signal to bypass softly.
* \param final_spectrum the spectrum to output from the plugin
* \param output_fft_buffer the unprocessed spectrum remaining in the fft buffer
* \param wet_dry mixing coefficient
* \param fft_size size of the fft
*/
void soft_bypass(Sprocessor *self)
{
    int k;

    for (k = 0; k < self->fft_size; k++)
    {
        self->processed_fft_spectrum[k] = (1.f - self->wet_dry) * self->original_fft_spectrum[k] + self->processed_fft_spectrum[k] * self->wet_dry;
    }
}

void sp_run(Sprocessor *self, float *fft_spectrum, float *enable)
{
    update_wetdry_target(self, enable);

    memcpy(self->original_fft_spectrum, fft_spectrum, sizeof(float) * self->fft_size);

    get_info_from_bins(self->power_spectrum, self->magnitude_spectrum,
                       self->phase_spectrum, self->fft_size_2,
                       self->fft_size, self->original_fft_spectrum);

    //Call the processing to be applied to the spectrum (in this case noise reduction)
    //This is temporary just to test   
    for (int i=0;i<self->fft_size;i++){
        self->processed_fft_spectrum[i] =  self->original_fft_spectrum[i] * 0.7f;
    }

    //If bypassed mix unprocessed and processed signal softly
    soft_bypass(self);

    //Copy the processed spectrum to fft_spectrum
    memcpy(fft_spectrum, self->processed_fft_spectrum, sizeof(float) * self->fft_size);
}