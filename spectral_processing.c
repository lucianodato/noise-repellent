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

#define SNR_INFLUENCE 1.0    //local SNR Influence for threshold scaing (from non linear subtraction)

//------------GAIN AND THRESHOLD CALCULATION---------------

//ADAPTIVE NOISE PROFILE
void spectral_gain_adaptive(float* fft_p2,
														float* fft_p2_prev_env,
												    float noise_thresholds_offset,
												    float* noise_thresholds_p2,
												    int fft_size_2,
														float release_coeff,
												    float* Gk){
	int k;
	float noise_thresholds_scaled[fft_size_2+1];

	//PREPROCESSING

	//Applying envelopes to signal power spectrum
	apply_envelope(fft_p2,
								 fft_p2_prev_env,
								 fft_size_2,
								 release_coeff);

	//OVERSUSTRACTION
	//Scale noise profile (equals applying an oversustraction factor in spectral sustraction)
	for (k = 0; k <= fft_size_2; k++) {
		noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset;
	}

	//GAIN CALCULATION
	power_subtraction(fft_size_2,
								    fft_p2,
								    noise_thresholds_scaled,
								    Gk);
}

//FOR MANUAL NOISE PROFILE
void spectral_gain_manual(float* fft_p2,
							    float* fft_p2_prev_tsmooth,
							    float* fft_p2_prev_env,
							    float* fft_p2_prev_tpres,
							    float time_smoothing,
									float artifact_control,
							    float noise_thresholds_offset,
							    float* noise_thresholds_p2,
							    int fft_size_2,
									float* prev_beta,
							    float* Gk,
									float release_coeff){

	int k;
	float noise_thresholds_scaled[fft_size_2+1];
	float Gk_wideband_gate;
	float original_spectrum[fft_size_2+1];
	float non_linear_snr;

	memcpy(original_spectrum,fft_p2,sizeof(float)*(fft_size_2+1));

	//PREPROCESSING

	//------SMOOTHING DETECTOR------

	//Applying envelopes to signal power spectrum
	apply_envelope(fft_p2,
								 fft_p2_prev_env,
								 fft_size_2,
								 release_coeff);

	/*Time smoothing between current and past power spectrum (similar effect to ephraim and malah)
		The best option here is to adaptively smooth 2D spectral components so it will require a biger buffer
		as suggested by Lukin in Suppression of Musical Noise Artifacts in Audio Noise Reduction by Adaptive 2D Filtering
	*/
	if (time_smoothing > 0.f){ //Issue 33 TODO
		spectrum_time_smoothing(fft_size_2,
														fft_p2_prev_tsmooth,
														fft_p2,
														time_smoothing);

		//Store previous power values for smoothing
		memcpy(fft_p2_prev_tsmooth,fft_p2,sizeof(float)*(fft_size_2+1));
	}

	//------OVERSUBTRACTION------

	//Scale noise thresholds (equals applying an oversubtraction factor in spectral subtraction)
	for (k = 0; k <= fft_size_2; k++) {
		//Adapting scaling of thresholds using local SNR as in Non linear subtraction
		//This could be adaptive using masking instead of local snr scaling TODO
		non_linear_snr = SNR_INFLUENCE + sqrtf(original_spectrum[k]/noise_thresholds_p2[k]);

		//Application of every scaling factor to noise thresholds
		noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset * non_linear_snr;
	}

	//------GAIN CALCULATION------

	//Spectral gate
	spectral_gating(fft_size_2,
									fft_p2,
									noise_thresholds_scaled,
									Gk);

	//------POSTPROCESSING GAINS------

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
}


//GAIN APPLICATION
void gain_application(int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
								      float* Gk,
											float whitening_factor,
											float tapering,
											float reduction_amount,
											float makeup_gain,
								      float wet_dry,
								      float noise_listen){

  int k;
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

	//POSTPROCESSING RESIDUAL
	//Whitening (residual spectrum more similar to white noise)
	//Tappering (preserves HF but reduces more lower ones)
	if(whitening_factor > 0.f) {
		whitening_and_tapering(residual_spectrum,whitening_factor,tapering,fft_size_2);
	}

	//OUTPUT RESULTS

  //Output processed signal or to noise only
  if (noise_listen == 0.f){
    //Output denoised result
		//Parametric application of the reduction
    for (k = 0; k <= fft_size_2; k++) {
      final_fft_buffer[k] =  denoised_fft_buffer[k] + residual_spectrum[k]*reduction_amount;
      if(k < fft_size_2)
        final_fft_buffer[fft_size-k] = denoised_fft_buffer[fft_size-k] + residual_spectrum[fft_size-k]*reduction_amount ;
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
    final_fft_buffer[k] *= makeup_gain;
    if(k < fft_size_2)
      final_fft_buffer[fft_size-k] *= makeup_gain;
  }

  //Smooth bypass
  for (k = 0; k <= fft_size_2; k++) {
    output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + final_fft_buffer[k] * wet_dry;
    if(k < fft_size_2)
      output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + final_fft_buffer[fft_size-k] * wet_dry;
  }
}
