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
} Spectral_Smoother;

//---------------TIME SMOOTHING--------------

/**
* Spectral smoothing proposed in 'Spectral subtraction with adaptive averaging of
* the gain function' but is not used yet.
* \param fft_size_2 half of the fft size
* \param spectrum the current power spectrum
* \param spectrum_prev the previous power spectrum
* \param noise_thresholds the noise thresholds estimated
* \param prev_beta beta corresponded to previos frame
* \param coeff reference smoothing value
*/
void spectrum_adaptive_time_smoothing(int fft_size_2, float *spectrum_prev, float *spectrum,
                                      float *noise_thresholds, float *prev_beta, float coeff)
{
    int k;
    float discrepancy, numerator = 0.f, denominator = 0.f;
    float beta_ts;
    float beta_smooth;
    float gamma_ts;

    for (k = 0; k <= fft_size_2; k++)
    {
        //These has to be magnitude spectrums
        numerator += fabs(spectrum[k] - noise_thresholds[k]);
        denominator += noise_thresholds[k];
    }
    //this is the discrepancy of the spectum
    discrepancy = numerator / denominator;
    //beta is the adaptive coefficient
    beta_ts = MIN(discrepancy, 1.f);

    //Gamma is the smoothing coefficient of the adaptive factor beta
    if (*prev_beta < beta_ts)
    {
        gamma_ts = 0.f;
    }
    else
    {
        gamma_ts = coeff;
    }

    //Smoothing beta
    beta_smooth = gamma_ts * *(prev_beta) + (1.f - gamma_ts) * beta_ts;

    //copy current value to previous
    *prev_beta = beta_smooth;

    //Apply the adaptive smoothed beta over the signal
    for (k = 0; k <= fft_size_2; k++)
    {
        spectrum[k] = (1.f - beta_smooth) * spectrum_prev[k] + beta_smooth * spectrum[k];
    }
}

void get_release_coeff(FFT_Denoiser *self, float release)
{
    //Parameters values
    /*exponential decay coefficients for envelopes and adaptive noise profiling
        These must take into account the hop size as explained in the following paper
        FFT-BASED DYNAMIC RANGE COMPRESSION*/
    if (release != 0.f) //This allows to turn off smoothing with 0 ms in order to use masking only
    {
        self->release_coeff = expf(-1000.f / (((release)*self->samp_rate) / self->hop));
    }
    else
    {
        self->release_coeff = 0.f; //This avoids incorrect results when moving sliders rapidly
    }
}

/**
* Spectral time smoothing by applying a release envelope. This seems to work better than * using time smoothing directly or McAulay & Malpass modification.
* \param spectrum the current power spectrum
* \param spectrum_prev the previous power spectrum
* \param N half of the fft size
* \param release_coeff release coefficient
*/
void apply_time_envelope(float *spectrum, float *spectrum_prev, float N, float release_coeff)
{
    int k;

    for (k = 0; k <= N; k++)
    {
        //It doesn't make much sense to have an attack slider when there is time smoothing
        if (spectrum[k] > spectrum_prev[k])
        {
            //Release (when signal is incrementing in amplitude)
            spectrum[k] = release_coeff * spectrum_prev[k] + (1.f - release_coeff) * spectrum[k];
        }
    }
}

memcpy(smoothed_spectrum, self->signal_spectrum, sizeof(float) * (self->half_fft_size + 1));

apply_time_envelope(smoothed_spectrum, smoothed_spectrum_prev, half_fft_size, release_coeff);

memcpy(smoothed_spectrum_prev, smoothed_spectrum, sizeof(float) * (half_fft_size + 1));