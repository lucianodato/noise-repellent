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
* \file fft_processor.c
* \author Luciano Dato
* \brief Contains an abstraction for a single fft spectrum processing
* \ It has an instance of fft denoiser wich contains everything to reduce
* \ the noise for current block fft
*/

//#include "fft_denoiser.c"
#include "extra_functions.c"

/**
* FFT processor struct.
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

    //Soft bypass
    float tau;            //time constant for soft bypass
    float wet_dry_target; //softbypass target for softbypass
    float wet_dry;        //softbypass coeff

    //FFT denoiser instance
    //FFTdenoiser *fft_denoiser;
} FFTprocessor;

/**
* Reset dynamic arrays to zero.
*/
void fft_p_reset(FFTprocessor *self)
{
    //Reset all arrays
    initialize_array(self->original_fft_spectrum, 0.f, self->fft_size);
    initialize_array(self->processed_fft_spectrum, 0.f, self->fft_size);
}

/**
* FFT processor initialization and configuration.
*/
FFTprocessor *
fft_p_init(int fft_size, int samp_rate, int hop)
{
    //Allocate object
    FFTprocessor *self = (FFTprocessor *)malloc(sizeof(FFTprocessor));

    //Configuration
    self->fft_size = fft_size;
    self->fft_size_2 = self->fft_size / 2;
    self->samp_rate = samp_rate;
    self->hop = hop;

    //spectrum allocation
    self->original_fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));
    self->processed_fft_spectrum = (float *)calloc((self->fft_size), sizeof(float));

    //soft bypass
    self->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->samp_rate));
    self->wet_dry = 0.f;

    //Initialize arrays with zeros
    fft_p_reset(self);

    //FFT denoiser initialization
    //fft_d_init(self->fft_denoiser);

    return self;
}

/**
* Free allocated memory.
*/
void fft_p_free(FFTprocessor *self)
{
    free(self->original_fft_spectrum);
    free(self->processed_fft_spectrum);
    //fft_d_free(self->fft_denoiser);
    free(self);
}

/**
* Updates the wet/dry mixing coefficient.
*/
void fft_p_update_wetdry_target(FFTprocessor *self, int enable)
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
void fft_p_soft_bypass(FFTprocessor *self)
{
    int k;

    for (k = 0; k < self->fft_size; k++)
    {
        self->processed_fft_spectrum[k] = (1.f - self->wet_dry) * self->original_fft_spectrum[k] + self->processed_fft_spectrum[k] * self->wet_dry;
    }
}

/**
* Runs the fft processing for current block.
*/
void fft_p_run(FFTprocessor *self, float *fft_spectrum, int enable)
{
    fft_p_update_wetdry_target(self, enable);

    memcpy(self->original_fft_spectrum, fft_spectrum, sizeof(float) * self->fft_size);

    /*Call the processing to be applied to the spectrum (in this case noise reduction)
      and store the processed spectrum in the processed fft spectrum
    */
    for (int i = 0; i < self->fft_size; i++)//This is just to test
    {
        self->processed_fft_spectrum[i] = self->original_fft_spectrum[i] * 0.7f;
    }
    //fft_d_process((self->fft_denoiser, self->original_fft_spectrum, adaptive_state, noise_learn_state));

    //If bypassed mix unprocessed and processed signal softly
    fft_p_soft_bypass(self);

    //Copy the processed spectrum to fft_spectrum
    memcpy(fft_spectrum, self->processed_fft_spectrum, sizeof(float) * self->fft_size);
}