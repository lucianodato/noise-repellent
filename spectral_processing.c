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

#include <float.h>
#include <math.h>

#include "estimate_noise_spectrum.c"
#include "denoise_gain.c"
#include "masking_thresholds.c"

void get_noise_thresholds(int auto_state,
                          int capture_state,
                          int fft_size_2,
                          float* power_spectrum,
                          float* magnitude_spectrum,
                          float* noise_thresholds,
                          float* auto_thresholds,
                          float* prev_noise_thresholds,
                          float* s_pow_spec,
                          float* prev_s_pow_spec,
                          float* p_min,
                          float* prev_p_min,
                          float* speech_p_p,
                          float* prev_speech_p_p,
                          float* window_count){

  //if slected auto estimate noise spectrum
  if(auto_state == 1.f) {
    auto_capture_noise(power_spectrum,//this is supposed to be the power spectrum in Loizou method
                       fft_size_2,
                       noise_thresholds,
                       auto_thresholds,
                       prev_noise_thresholds,
                       s_pow_spec,
                       prev_s_pow_spec,
                       p_min,
                       prev_p_min,
                       speech_p_p,
                       prev_speech_p_p);
  }

  //If selected estimate noise spectrum based on selected portion of signal
  if(capture_state == 1.f) { //MANUAL
    get_noise_statistics(magnitude_spectrum,
                         fft_size_2,
                         noise_thresholds,
                         window_count);
  }
}

//------------GAIN AND THRESHOLD CALCULATION---------------

void spectral_gain_computing(float* bark_z,
                    float* fft_p2,
                    float* fft_p2_prev,
                    float* fft_p2_smooth,
                    float* fft_magnitude,
                    float* fft_magnitude_prev,
                    float* fft_magnitude_smooth,
                    float time_smoothing,
                    float* noise_thresholds,
                    int fft_size_2,
                    float* alpha,
                    //float* beta,
                    float max_masked,
                    float min_masked,
                    float reduction_strenght,
                    float* Gk,
                    float* Gk_prev,
                    float masking,
                    float frequency_smoothing){

  //PREPROCESSING

  //CALCULATION OF ALPHA AND BETA WITH MASKING THRESHOLDS USING VIRAG METHOD
  if(masking == 1.f){
    compute_alpha_and_beta(bark_z,
                           fft_p2,
                           fft_magnitude,
                           noise_thresholds,
                           fft_size_2,
                           alpha,
                           //beta,
                           max_masked,
                           min_masked);
  }

  //SMOOTHING
  //Time smoothing between current and past power spectrum and magnitude spectrum
  spectrum_exponential_smoothing(fft_size_2,
                                 fft_p2_prev,
                                 fft_p2,
                                 fft_p2_smooth,
                                 time_smoothing);

  spectrum_exponential_smoothing(fft_size_2,
                                 fft_magnitude_prev,
                                 fft_magnitude,
                                 fft_magnitude_smooth,
                                 time_smoothing);

  //GAIN CALCULATION
  if(masking == 1.f){
    //Parametric Generalized Spectral Sustraction
    //(Power Sustraction with variable alpha)
    denoise_gain_gss_fixed_beta(reduction_strenght,
                                fft_size_2,
                                alpha,
                                fft_magnitude_smooth,
                                noise_thresholds,
                                Gk,
                                Gk_prev);
  } else {
    //Power Sustraction
    denoise_gain_gss(reduction_strenght,
                     fft_size_2,
                     ALPHA_GSS,//alpha
                     BETA_GSS,//beta
                     fft_magnitude_smooth,
                     noise_thresholds,
                     Gk);
  }

  //FREQUENCY SMOOTHING OF GAINS
  spectral_smoothing_MA(Gk,
                        frequency_smoothing,
                        fft_size_2);
}

//GAIN APPLICATION
void gain_application(float amount_of_reduction,
                      int fft_size_2,
                      int fft_size,
                      float* output_fft_buffer,
                      float* Gk,
                      float wet_dry,
                      float residual_whitening,
                      float noise_listen){

  int k;
  float reduction_coeff = (1.f/from_dB(amount_of_reduction));
  float denoised_fft_buffer[fft_size];
  float residual_spectrum[fft_size];
  float tappering_filter[fft_size_2+1];

  //Apply the computed gain to the signal and store it in denoised array
  for (k = 0; k <= fft_size_2; k++) {
    denoised_fft_buffer[k] = output_fft_buffer[k] * Gk[k];
    if(k < fft_size_2)
      denoised_fft_buffer[fft_size-k] = output_fft_buffer[fft_size-k] * Gk[k];
  }

  //Residual signal
  for (k = 0; k <= fft_size_2; k++) {
   residual_spectrum[k] = output_fft_buffer[k] - denoised_fft_buffer[k];
   if(k < fft_size_2)
    residual_spectrum[fft_size-k] = output_fft_buffer[fft_size-k] - denoised_fft_buffer[fft_size-k];
  }

  //Residual signal Whitening and tappering
  if(residual_whitening != 0.f) {
    whitening_of_spectrum(residual_spectrum,residual_whitening,fft_size_2);
    tappering_filter_calc(tappering_filter,(fft_size_2+1),residual_whitening);
    apply_tappering_filter(residual_spectrum,tappering_filter,fft_size_2);
  }

  //Listen to cleaned signal or to noise only
  if (noise_listen == 0.f){
    //Mix residual and processed (Parametric way of noise reduction)
    for (k = 0; k <= fft_size_2; k++) {
      output_fft_buffer[k] =  (1.f-wet_dry) * output_fft_buffer[k] + (denoised_fft_buffer[k] + residual_spectrum[k]*reduction_coeff) * wet_dry;
      if(k < fft_size_2)
        output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + (denoised_fft_buffer[fft_size-k] + residual_spectrum[fft_size-k]*reduction_coeff) * wet_dry;
    }
  } else {
    //Output noise only
    for (k = 0; k <= fft_size_2; k++) {
      output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + residual_spectrum[k] * wet_dry;
      if(k < fft_size_2)
        output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + residual_spectrum[fft_size-k] * wet_dry;
    }
  }
}
