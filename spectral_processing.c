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
														 float artifact_control,
												     float noise_thresholds_offset,
														 float auto_state,
												     float* noise_thresholds_p2,
												     int fft_size_2,
												     int fft_size,
												     float* Gk,
													 	 float* Gk_prev,
													 	 float* Gk_prev_wide,
														 float attack_coeff,
														 float release_coeff,
													 	 float knee_width){

	//PREPROCESSING
	int k;
	float noise_thresholds_scaled[fft_size_2+1];
	//float Gk_wideband_gate[fft_size_2+1];
	float Gk_power_sustraction[fft_size_2+1];
	float Gk_spectral_gate[fft_size_2+1];

	//Scale noise profile (equals applying an oversustraction factor in spectral sustraction)
	//This must be adaptive using masking or local snr strategy
	// if (auto_state != 1.f){
		for (k = 0; k <= fft_size_2; k++) {
			noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset * (1.f + sqrtf(fft_p2[k]/noise_thresholds_p2[k]));
		}
	// }else{
	// 	//Prevent local snr scaling when using adaptive profiling
	// 	for (k = 0; k <= fft_size_2; k++) {
	// 		noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset;
	// 	}
	// }


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
							     Gk_power_sustraction);


	spectral_gating(fft_size_2,
		    attack_coeff,
				release_coeff,
				knee_width,
		    fft_p2,
		    noise_thresholds_scaled,
		    Gk_spectral_gate,
		    Gk_prev);

	// wideband_gating(fft_size_2,
	// 	    attack_coeff,
	// 			release_coeff,
	// 			knee_width,
	// 	    fft_p2,
	// 	    noise_thresholds_scaled,
	// 	    Gk_wideband_gate,
	// 	    Gk_prev_wide);

	//Artifact control (interpolation between power sustraction and gating strategies)
	for (k = 0; k <= fft_size_2; k++) {
		Gk[k] = (1.f - artifact_control)*Gk_power_sustraction[k] + artifact_control*Gk_spectral_gate[k];
		//Gk[k] = (1.f - artifact_control)*Gk_spectral_gate[k] + artifact_control*Gk_wideband_gate[k];
	}
}

//GAIN APPLICATION
void gain_application(float amount_of_reduction,
								      int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
								      float* Gk,
											float whitening_factor,
											float tapering,
								      float makeup_gain,
								      float wet_dry,
								      float noise_listen){

  int k;
  float reduction_amount = from_dB(-1.f*amount_of_reduction);
	float gain = from_dB(makeup_gain);
  float residual_spectrum[fft_size];
  float tapering_filter[fft_size];
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

	//Whitening and tappering
	if(whitening_factor > 0.f) {
		whitening_of_spectrum(residual_spectrum,whitening_factor,fft_size_2);
		if(tapering > 0.f){
			tapering_filter_calc(tapering_filter,(fft_size_2+1));
			apply_tapering_filter(residual_spectrum,tapering_filter,fft_size_2);
		}
	}

  //Listen to processed signal or to noise only
  if (noise_listen == 0.f){
    //Mix residual and processed (Parametric way of noise reduction)
    for (k = 0; k <= fft_size_2; k++) {
      final_fft_buffer[k] =  denoised_fft_buffer[k] + (residual_spectrum[k]*reduction_amount);
      if(k < fft_size_2)
        final_fft_buffer[fft_size-k] = denoised_fft_buffer[fft_size-k] + (residual_spectrum[fft_size-k]*reduction_amount);
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
    final_fft_buffer[k] *= gain;
    if(k < fft_size_2)
      final_fft_buffer[fft_size-k] *= gain;
  }

  //Smooth bypass
  for (k = 0; k <= fft_size_2; k++) {
    output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + final_fft_buffer[k] * wet_dry;
    if(k < fft_size_2)
      output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + final_fft_buffer[fft_size-k] * wet_dry;
  }
}
