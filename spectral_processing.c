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

//------------GAIN AND THRESHOLD CALCULATION---------------

void spectral_gain_computing(float* fft_p2,
                             float* fft_p2_prev,
                             float* fft_magnitude,
                             float* fft_magnitude_prev,
                             float time_smoothing,
                             float reduction_scale,
                             float* noise_thresholds_p2,
                             float* noise_thresholds_magnitude,
                             int fft_size_2,
                             float reduction_strenght,
                             float* Gk,
                             float frequency_smoothing){

  //PREPROCESSING

  //SMOOTHING
  //Time smoothing between current and past power spectrum and magnitude spectrum
  if (time_smoothing > 0.f){
    spectrum_time_smoothing(fft_size_2,
                            fft_p2_prev,
                            fft_p2,
                            time_smoothing);

    spectrum_time_smoothing(fft_size_2,
                            fft_magnitude_prev,
                            fft_magnitude,
                            time_smoothing);
  }

  //GAIN CALCULATION
  //Non linear Power Sustraction
  nonlinear_power_sustraction(reduction_strenght,
                   reduction_scale,
                   fft_size_2,
                   fft_magnitude,
                   noise_thresholds_magnitude,
                   Gk);

  //FREQUENCY SMOOTHING OF GAINS
  if (frequency_smoothing > 0.f){
    spectral_smoothing_MA(Gk,
                          frequency_smoothing,
                          fft_size_2);
  }
}

//GAIN APPLICATION
void gain_application(float amount_of_reduction,
                      int fft_size_2,
                      int fft_size,
                      float* output_fft_buffer,
                      float* Gk,
                      float makeup_gain,
                      float wet_dry,
                      float residual_whitening,
                      float noise_listen){

  int k;
  float reduction_coeff = from_dB(-1.f*sanitize_denormal(amount_of_reduction));
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
  if(residual_whitening > 0.f) {
    whitening_of_spectrum(residual_spectrum,residual_whitening,fft_size_2);
    tappering_filter_calc(tappering_filter,(fft_size_2+1));
    apply_tappering_filter(residual_spectrum,tappering_filter,fft_size_2);
  }

  //Listen to cleaned signal or to noise only
  if (noise_listen == 0.f){
    //Mix residual and processed (Parametric way of noise reduction)
    for (k = 0; k <= fft_size_2; k++) {
      output_fft_buffer[k] =  (1.f-wet_dry) * output_fft_buffer[k] + from_dB(makeup_gain) * (denoised_fft_buffer[k] + residual_spectrum[k]*reduction_coeff) * wet_dry;
      if(k < fft_size_2)
        output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + from_dB(makeup_gain) * (denoised_fft_buffer[fft_size-k] + residual_spectrum[fft_size-k]*reduction_coeff) * wet_dry;
    }
  } else {
    //Output noise only
    for (k = 0; k <= fft_size_2; k++) {
      output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + from_dB(makeup_gain) * residual_spectrum[k] * wet_dry;
      if(k < fft_size_2)
        output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + from_dB(makeup_gain) * residual_spectrum[fft_size-k] * wet_dry;
    }
  }
}
