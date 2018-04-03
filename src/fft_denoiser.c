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
* \file fft_denoiser.c
* \author Luciano Dato
* \brief Contains a the spectral noise reduction abstraction
*/

#include "gain_estimator.c"
#include "noise_estimator.c"

/**
* FFT denoising struct.
*/
typedef struct
{
    //General parameters
    int fft_size;
    int fft_size_2;
    int samp_rate;
    int hop;

    //Reduction gains
    float *Gk; //definitive gain

    //Ensemble related
    //Spectrum
    float *power_spectrum;
    float *phase_spectrum;
    float *magnitude_spectrum;

    float *residual_spectrum;
    float *denoised_spectrum;
    float *final_spectrum;

    Gestimator *gain_estimation;
    Nestimator *noise_estimation;

    //Parameters for the algorithm (user input)
    float amount_of_reduction;     //Amount of noise to reduce in dB

    //Algorithm exta variables
    float reduction_coeff;            //Gain to apply to the residual noise
    float amount_of_reduction_linear; //Reduction amount linear value
    float whitening_factor;           //Whitening amount of the reduction

    //whitening related
    float *residual_max_spectrum;
    float max_decay_rate;
    float whitening_window_count;
} FFTdenoiser;

/**
* To reset the noise profile and set every value to default one.
*/
static void
fft_d_reset_noise_profile(FFTdenoiser *self)
{
    initialize_array(self->Gk, 1.f, self->fft_size);

    initialize_array(self->power_spectrum, 0.f, self->fft_size_2 + 1);
    initialize_array(self->magnitude_spectrum, 0.f, self->fft_size_2 + 1);
    initialize_array(self->phase_spectrum, 0.f, self->fft_size_2 + 1);

    initialize_array(self->residual_max_spectrum, 0.f, self->fft_size);
    self->whitening_window_count = 0.f;
}

void fft_d_init(FFTdenoiser *self, int fft_size, int samp_rate, int hop, float amount_of_reduction,
                float whitening_factor_pc)
{
    //Configuration
    self->fft_size = fft_size;
    self->fft_size_2 = self->fft_size / 2;
    self->samp_rate = samp_rate;
    self->hop = hop;
    self->amount_of_reduction_linear = from_dB(-1.f * *(amount_of_reduction));
    self->whitening_factor = *(whitening_factor_pc) / 100.f;

    //reduction gains related
    self->Gk = (float *)calloc((self->fft_size), sizeof(float));

    //whitening related
    self->residual_max_spectrum = (float *)calloc((self->fft_size), sizeof(float));
    self->max_decay_rate = expf(-1000.f / (((WHITENING_DECAY_RATE)*samp_rate) / self->hop));
    self->whitening_window_count = 0.f;

    //final ensemble related
    self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
    self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));
    self->final_spectrum = (float *)calloc((self->fft_size), sizeof(float));

    //Arrays for getting bins info
    self->power_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));
    self->magnitude_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));
    self->phase_spectrum = (float *)malloc((self->fft_size_2 + 1) * sizeof(float));

    //Set initial gain as unity for the positive part
    initialize_array(self->Gk, 1.f, self->fft_size);
}

void fft_p_free(FFTdenoiser *self)
{
    free(self->power_spectrum);
    free(self->magnitude_spectrum);
    free(self->phase_spectrum);
    free(self);
}

void fft_d_process(FFTdenoiser *self, float* original_fft_spectrum, int adaptive_state,
                   int noise_learn_state)
{
    //First get the power, magnitude and phase spectrum 
    get_info_from_bins(self->power_spectrum, self->magnitude_spectrum,
                       self->phase_spectrum, self->fft_size_2,
                       self->fft_size, self->original_fft_spectrum);

    //If the spectrum is not silence
    if (!is_empty(self->power_spectrum, self->fft_size_2))
    {
        //If adaptive noise is selected the noise is adapted in time
        if (*(adaptive_state) == 1.f)
        {
            //This has to be revised(issue 8 on github)
            adapt_noise(self->power_spectrum, self->fft_size_2, self->noise_thresholds_p2,
                        self->auto_thresholds, self->prev_noise_thresholds,
                        self->s_pow_spec, self->prev_s_pow_spec, self->p_min,
                        self->prev_p_min, self->speech_p_p, self->prev_speech_p_p);

            self->noise_thresholds_availables = true;
        }

        /*If selected estimate noise spectrum is based on selected portion of signal
            *do not process the signal
            */
        if (*(noise_learn_state) == 1.f)
        { //MANUAL

            //Increase window count for rolling mean
            self->noise_window_count++;

            get_noise_statistics(self->power_spectrum, self->fft_size_2,
                                 self->noise_thresholds_p2, self->noise_window_count);

            self->noise_thresholds_availables = true;
        }
        else
        {
            //If there is a noise profile reduce noise
            if (self->noise_thresholds_availables == true)
            {
                //Detector smoothing and oversubtraction
                preprocessing(self->thresholds_offset_linear, self->power_spectrum,
                              self->noise_thresholds_p2, self->noise_thresholds_scaled,
                              self->smoothed_spectrum, self->smoothed_spectrum_prev,
                              self->fft_size_2, self->bark_z, self->absolute_thresholds,
                              self->SSF, self->release_coeff,
                              self->spreaded_unity_gain_bark_spectrum,
                              self->spl_reference_values, self->alpha_masking,
                              self->beta_masking, *(self->masking), *(self->adaptive_state),
                              self->amount_of_reduction_linear, self->transient_preserv_prev,
                              &self->tp_window_count, &self->tp_r_mean,
                              &self->transient_present, *(self->transient_protection));

                //Supression rule
                spectral_gain(self->power_spectrum, self->noise_thresholds_p2,
                              self->noise_thresholds_scaled, self->smoothed_spectrum,
                              self->fft_size_2, *(self->adaptive_state), self->Gk,
                              *(self->transient_protection), self->transient_present);

                //apply gains
                denoised_calulation(self->fft_size, self->fft_spectrum,
                                    self->denoised_spectrum, self->Gk);

                //residual signal
                residual_calulation(self->fft_size, self->fft_spectrum,
                                    self->residual_spectrum, self->denoised_spectrum,
                                    self->whitening_factor, self->residual_max_spectrum,
                                    &self->whitening_window_count, self->max_decay_rate);

                //Ensemble the final spectrum using residual and denoised
                final_spectrum_ensemble(self->fft_size, self->final_spectrum,
                                        self->residual_spectrum,
                                        self->denoised_spectrum,
                                        self->amount_of_reduction_linear,
                                        *(self->residual_listen));
            }
        }
    }
}

/**
* Applies the filter to the complex spectrum and gets the clean signal.
* \param fft_size size of the fft
* \param output_fft_buffer the unprocessed spectrum remaining in the fft buffer
* \param denoised_spectrum the spectrum of the cleaned signal
* \param Gk is the filter computed by the supression rule for each bin of the spectrum
*/
void denoised_calulation(int fft_size, float *output_fft_buffer,
						 float *denoised_spectrum, float *Gk)
{
	int k;

	//Apply the computed gain to the signal and store it in denoised array
	for (k = 0; k < fft_size; k++)
	{
		denoised_spectrum[k] = output_fft_buffer[k] * Gk[k];
	}
}

/**
* Gets the residual signal of the reduction.
* \param fft_size size of the fft
* \param output_fft_buffer the unprocessed spectrum remaining in the fft buffer
* \param denoised_spectrum the spectrum of the cleaned signal
* \param whitening_factor the mix coefficient between whitened and not whitened residual spectrum
* \param residual_max_spectrum contains the maximun temporal value in each residual bin
* \param whitening_window_count counts frames to distinguish the first from the others
* \param max_decay_rate coefficient that sets the memory for each temporal maximun
*/
void residual_calulation(int fft_size, float *output_fft_buffer,
						 float *residual_spectrum, float *denoised_spectrum,
						 float whitening_factor, float *residual_max_spectrum,
						 float *whitening_window_count, float max_decay_rate)
{

	int k;

	//Residual signal
	for (k = 0; k < fft_size; k++)
	{
		residual_spectrum[k] = output_fft_buffer[k] - denoised_spectrum[k];
	}

	////////////POSTPROCESSING RESIDUAL
	//Whitening (residual spectrum more similar to white noise)
	if (whitening_factor > 0.f)
	{
		spectral_whitening(residual_spectrum, whitening_factor, fft_size,
						   residual_max_spectrum, whitening_window_count, max_decay_rate);
	}
	////////////
}

/**
* Mixes the cleaned signal with the residual taking into account the reduction configured
* by the user. Outputs the final signal or the residual only.
* \param fft_size size of the fft
* \param final_spectrum the spectrum to output from the plugin
* \param residual_spectrum the spectrum of the reduction residual
* \param denoised_spectrum the spectrum of the cleaned signal
* \param reduction_amount the amount of dB power to reduce setted by the user
* \param noise_listen control variable that decides whether to output the mixed noise reduced signal or the residual only
*/
void final_spectrum_ensemble(int fft_size, float *final_spectrum,
							 float *residual_spectrum, float *denoised_spectrum,
							 float reduction_amount, float noise_listen)
{
	int k;

	//OUTPUT RESULTS using smooth bypass and parametric subtraction
	if (noise_listen == 0.f)
	{
		//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k < fft_size; k++)
		{
			final_spectrum[k] = denoised_spectrum[k] + residual_spectrum[k] * reduction_amount;
		}
	}
	else
	{
		//Output noise only
		for (k = 0; k < fft_size; k++)
		{
			final_spectrum[k] = residual_spectrum[k];
		}
	}
}