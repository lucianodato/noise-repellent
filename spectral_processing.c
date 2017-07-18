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

#define SNR_INFLUENCE 1.0f    //local SNR Influence for threshold scaing (from non linear sustraction)

//------------GAIN AND THRESHOLD CALCULATION---------------

//MANUAL NOISE PROFILE
void spectral_gain_manual(float* fft_p2,
											    float* fft_p2_prev_tsmooth,
											    float* fft_p2_prev_env,
											    float* fft_p2_prev_tpres,
													float amount_of_reduction,
		 											float whitening_factor,
		 											float tapering,
											    float time_smoothing,
													float artifact_control,
											    float noise_thresholds_offset,
											    float transient_preservation_switch,
											    float* noise_thresholds_p2,
											    int fft_size_2,
											    float* Gk,
												 	float* Gk_prev,
												 	float* Gk_prev_wide,
													float release_coeff){

	int k;
	float noise_thresholds_scaled[fft_size_2+1];
	float Gk_wideband_gate;
	float original_spectrum[fft_size_2+1];
	memcpy(original_spectrum,fft_p2,sizeof(float)*(fft_size_2+1));
	float transient_preservation_coeff = 1.f;
	float oversustraction_factor;
	float reduction_amount = from_dB(-1.f*amount_of_reduction);
	//float tapering_filter[fft_size_2+1];

	//PREPROCESSING

	//SMOOTHING DETECTOR

	//Applying envelopes to signal power spectrum
	apply_envelope(fft_p2,
								 fft_p2_prev_env,
								 fft_size_2,
								 release_coeff);

	//Time smoothing between current and past power spectrum (similar effect to ephraim and malah)
	if (time_smoothing > 0.f){
		spectrum_time_smoothing(fft_size_2,
														fft_p2_prev_tsmooth,
														fft_p2,
														time_smoothing);

		//Store previous power values for smoothing
		memcpy(fft_p2_prev_tsmooth,fft_p2,sizeof(float)*(fft_size_2+1));
	}

	//OVERSUSTRACTION

	//Transient preservation using onset detection
	if(transient_preservation_switch > 0.f){
		transient_preservation_coeff = transient_preservation(original_spectrum,
																													fft_p2_prev_tpres,
																													fft_size_2);

		//printf("%f\n", spectral_flux_value );

		memcpy(fft_p2_prev_tpres,original_spectrum,sizeof(float)*(fft_size_2+1));
	}

	//Scale noise profile (equals applying an oversustraction factor in spectral sustraction)
	for (k = 0; k <= fft_size_2; k++) {
		//This could be adaptive using masking instead of local snr scaling TODO
		oversustraction_factor = noise_thresholds_offset * transient_preservation_coeff * (SNR_INFLUENCE + sqrtf(fft_p2[k]/noise_thresholds_p2[k]));
		noise_thresholds_scaled[k] = noise_thresholds_p2[k] * oversustraction_factor;
	}

	//GAIN CALCULATION

	//POWER SUSTRACTION FILTER
	// power_sustraction(fft_size_2,
	// 							    fft_p2,
	// 									reduction_amount,
	// 							    noise_thresholds_scaled,
	// 							    Gk_power_sustraction);
	//
	// memcpy(Gk,Gk_power_sustraction,sizeof(float)*(fft_size_2+1));

	//SPECTRAL GATE
	spectral_gating(fft_size_2,
									fft_p2,
									reduction_amount,
									noise_thresholds_scaled,
									Gk);

	//POSTPROCESSING GAINS

	//Artifact control (applying wideband gating between pauses)
	if(artifact_control > 0.f){
		wideband_gating(fft_size_2,
										fft_p2,
										noise_thresholds_scaled,
										&Gk_wideband_gate);

		for (k = 0; k <= fft_size_2; k++) {
			//Only apply wideband gate when signal is below the threshold and scale is set by the user
			if (Gk_wideband_gate < 1.f)
				Gk[k] = (1.f-artifact_control)*Gk[k] +  artifact_control*Gk_wideband_gate;
		}
	}

	// //Whitening and tappering POSTPROCESSING
	// if(whitening_factor > 0.f) {
	// 	whitening_of_spectrum(residual_spectrum,whitening_factor,fft_size_2);
	// 	if(tapering > 0.f){
	// 		tapering_filter_calc(tapering_filter,(fft_size_2+1));
	// 		apply_tapering_filter(residual_spectrum,tapering_filter,fft_size_2);
	// 	}
	// }

}

//ADAPTIVE NOISE PROFILE
void spectral_gain_adaptive(float* fft_p2,
												    float noise_thresholds_offset,
												    float* noise_thresholds_p2,
												    int fft_size_2,
														float amount_of_reduction,
												    float* Gk){
	int k;
	float noise_thresholds_scaled[fft_size_2+1];
	float reduction_amount = from_dB(-1.f*amount_of_reduction);

	//PREPROCESSING

	//OVERSUSTRACTION
	//Scale noise profile (equals applying an oversustraction factor in spectral sustraction)
	for (k = 0; k <= fft_size_2; k++) {
		noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset;
	}

	//GAIN CALCULATION
	power_sustraction(fft_size_2,
								    fft_p2,
										reduction_amount,
								    noise_thresholds_scaled,
								    Gk);
}

//GAIN APPLICATION
void gain_application(int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
								      float* Gk,
											float makeup_gain,
								      float wet_dry,
								      float noise_listen){

  int k;
  float residual_spectrum[fft_size];
  float denoised_fft_buffer[fft_size];
  float final_fft_buffer[fft_size];
	float gain = from_dB(makeup_gain);

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

  //Listen to processed signal or to noise only
  if (noise_listen == 0.f){
    //Mix residual and processed (Parametric way of noise reduction)
    for (k = 0; k <= fft_size_2; k++) {
      final_fft_buffer[k] =  denoised_fft_buffer[k];
      if(k < fft_size_2)
        final_fft_buffer[fft_size-k] = denoised_fft_buffer[fft_size-k];
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
