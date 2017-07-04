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
												     float time_smoothing,
												     float noise_thresholds_offset,
												     float* noise_thresholds_p2,
												     int fft_size_2,
												     int fft_size,
												     float* Gk){

	//PREPROCESSING
	int k;
	float noise_thresholds_scaled[fft_size_2+1];

	//Scale noise profile (equals applying an oversustraction factor)
	for (k = 0; k <= fft_size_2; k++) {
		noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset;
	}

	//SMOOTHING
	//Time smoothing between current and past power spectrum and magnitude spectrum
	if (time_smoothing > 0.f){
		spectrum_time_smoothing(fft_size_2,
												    fft_p2_prev,
												    fft_p2,
												    time_smoothing);
	}

	//GAIN CALCULATION
	power_sustraction(fft_size_2,
							     fft_p2,
							     noise_thresholds_scaled,
							     Gk);
}

//GAIN APPLICATION
void gain_application(float amount_of_reduction,
								      int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
								      float* Gk,
											float* whitening_influence,
											float whitening_factor,
								      float makeup_gain,
								      float wet_dry,
								      float noise_listen){

  int k;
  float reduction_coeff = from_dB(-1.f*amount_of_reduction);
	float freq_scaling;
	float reduction_influence[fft_size];
  float residual_spectrum[fft_size];
  float denoised_fft_buffer[fft_size];
  float final_fft_buffer[fft_size];

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

	//Apply scaling to the reduction amount in such way that the residual is more like white noise
	for (k = 0; k <= fft_size_2; k++) {
		freq_scaling = whitening_influence[k]*reduction_coeff;//This isn't quite correct
		//This interpolates between whitening and no whitening reduction amount
		reduction_influence[k] =  reduction_coeff*(1.f-whitening_factor) + whitening_factor*freq_scaling;
		if(k < fft_size_2)
			reduction_influence[fft_size-k] = reduction_influence[k];//mirroring frequencies
	}

  //Listen to processed signal or to noise only
  if (noise_listen == 0.f){
    //Mix residual and processed (Parametric way of noise reduction)
    for (k = 0; k <= fft_size_2; k++) {
      final_fft_buffer[k] =  denoised_fft_buffer[k] + residual_spectrum[k]*reduction_influence[k];
      if(k < fft_size_2)
        final_fft_buffer[fft_size-k] = denoised_fft_buffer[fft_size-k] + residual_spectrum[fft_size-k]*reduction_influence[fft_size-k];
    }
  } else {
    //Output noise only
    for (k = 0; k <= fft_size_2; k++) {
      final_fft_buffer[k] = residual_spectrum[k];
      if(k < fft_size_2)
        final_fft_buffer[fft_size-k] = residual_spectrum[fft_size-k];
    }
  }

  //Applying make up gain
  for (k = 0; k <= fft_size_2; k++) {
    final_fft_buffer[k] *= from_dB(makeup_gain);
    if(k < fft_size_2)
      final_fft_buffer[fft_size-k] *= from_dB(makeup_gain);
  }

  //Smooth bypass
  for (k = 0; k <= fft_size_2; k++) {
    output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + final_fft_buffer[k] * wet_dry;
    if(k < fft_size_2)
      output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + final_fft_buffer[fft_size-k] * wet_dry;
  }
}
