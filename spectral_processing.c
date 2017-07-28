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

//#define SNR_INFLUENCE 1.0    //local SNR Influence for threshold scaing (from non linear subtraction)

#define MIX_COEFF 0.2		//wideband gate MIX_COEFF wide<->spectral

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

	memcpy(fft_p2_prev_env,fft_p2,sizeof(float)*(fft_size_2+1));


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
							    float* fft_p2_prev_env,
							    float noise_thresholds_offset,
							    float* noise_thresholds_p2,
							    int fft_size_2,
							    float* Gk,
									float release_coeff){

	int k;
	float noise_thresholds_scaled[fft_size_2+1];

	//PREPROCESSING

	//------OVERSUBTRACTION------

	//Scale noise thresholds (equals applying an oversubtraction factor in spectral subtraction)
	for (k = 0; k <= fft_size_2; k++) {
		//Application of every scaling factor to noise thresholds
		noise_thresholds_scaled[k] = noise_thresholds_p2[k] * noise_thresholds_offset;
	}

	//------SMOOTHING DETECTOR------

	/*Time smoothing between current and past Gk (similar effect to ephraim and malah)
		Here is done by applying a release envelope to signal power spectrum
		The best option here is to adaptively smooth 2D spectral components so it will require a biger buffer
		as suggested by Lukin in Suppression of Musical Noise Artifacts in Audio Noise Reduction by Adaptive 2D Filtering
	*/
	apply_envelope(fft_p2,
								 fft_p2_prev_env,
								 fft_size_2,
								 release_coeff);

	memcpy(fft_p2_prev_env,fft_p2,sizeof(float)*(fft_size_2+1));

	//------GAIN CALCULATION------

	//Spectral gate
	spectral_gating(fft_size_2,
									fft_p2,
									noise_thresholds_scaled,
									Gk);

}


//GAIN APPLICATION
void gain_application(int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
											float* input_fft_buffer_ps,
											float* input_fft_buffer_g,
											float* output_fft_buffer_ps,
											float* output_fft_buffer_g,
											fftwf_plan* forward_g,
											fftwf_plan* backward_g,
											fftwf_plan* forward_ps,
											fftwf_plan* backward_ps,
								      float* Gk,
											float whitening_factor,
											float tapering,
											float reduction_amount,
											float makeup_gain,
								      float wet_dry,
								      float noise_listen){

  int k;
  float residual_spectrum[fft_size];
  float denoised_spectrum[fft_size];
  float postfilter[fft_size];
  float Gk_mirrored[fft_size];

	//mirrored gain array to use a power of 2 fft transform
	for (k = 0; k <= fft_size_2; k++) {
		Gk_mirrored[k] = Gk[k];
		if(k < fft_size_2)
			Gk_mirrored[fft_size-k] = Gk[k];
	}

  ////////////POSTPROCESSING GAIN SMOOTHING USING A POSTFILTER
	//Compute the filter
	compute_post_filter(fft_size_2,
										fft_size,
										output_fft_buffer,
										0.4f,//snr_threshold,
										2.f,//amount_smoothing,
										postfilter,
										Gk_mirrored);

	//Convolution using fft transform

	//Copy to fft buffers
	memcpy(input_fft_buffer_ps,postfilter,fft_size*sizeof(float));
	memcpy(input_fft_buffer_g,Gk_mirrored,fft_size*sizeof(float));

	//FFT Analysis
	fftwf_execute(*forward_ps);
	fftwf_execute(*forward_g);

	//Multiply with the filter computed
	for (k = 0; k <= fft_size_2; k++) {
    output_fft_buffer_g[k] *= output_fft_buffer_ps[k];
    if(k < fft_size_2)
      output_fft_buffer_g[fft_size-k] *= output_fft_buffer_ps[fft_size-k];
  }

	//FFT Synthesis
	fftwf_execute(*backward_ps);
	fftwf_execute(*backward_g);

	//Normalizing
	for (k = 0; k < fft_size; k++){
		input_fft_buffer_ps[k] = input_fft_buffer_ps[k] / fft_size;
		input_fft_buffer_g[k] = input_fft_buffer_g[k] / fft_size;
	}

	//Copy to orginal arrays
	memcpy(Gk_mirrored,input_fft_buffer_g,fft_size*sizeof(float));
	memcpy(postfilter,input_fft_buffer_ps,fft_size*sizeof(float));
	////////////

  //Apply the computed gain to the signal and store it in denoised array
  for (k = 0; k <= fft_size_2; k++) {
    denoised_spectrum[k] = output_fft_buffer[k] * Gk_mirrored[k];
    if(k < fft_size_2)
      denoised_spectrum[fft_size-k] = output_fft_buffer[fft_size-k] * Gk_mirrored[fft_size-k];
  }

  //Residual signal
  for (k = 0; k <= fft_size_2; k++) {
   residual_spectrum[k] = output_fft_buffer[k] - denoised_spectrum[k];
   if(k < fft_size_2)
    residual_spectrum[fft_size-k] = output_fft_buffer[fft_size-k] - denoised_spectrum[fft_size-k];
  }

	////////////POSTPROCESSING RESIDUAL
	//Whitening (residual spectrum more similar to white noise)
	//Tappering (preserves HF but reduces more lower ones)
	if(whitening_factor > 0.f) {
		whitening_and_tapering(residual_spectrum,whitening_factor,tapering,fft_size_2);
	}
	////////////

	//OUTPUT RESULTS using smooth bypass and parametric sustraction
	if (noise_listen == 0.f){
	//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k <= fft_size_2; k++) {
			output_fft_buffer[k] =  (1.f-wet_dry) * output_fft_buffer[k] + makeup_gain * (denoised_spectrum[k] + residual_spectrum[k]*reduction_amount) * wet_dry;
			if(k < fft_size_2)
				output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + makeup_gain * (denoised_spectrum[fft_size-k] + residual_spectrum[fft_size-k]*reduction_amount) * wet_dry;
		}
	} else {
		//Output noise only
		for (k = 0; k <= fft_size_2; k++) {
			output_fft_buffer[k] = (1.f-wet_dry) * output_fft_buffer[k] + makeup_gain * residual_spectrum[k] * wet_dry;
			if(k < fft_size_2)
				output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + makeup_gain * residual_spectrum[fft_size-k] * wet_dry;
		}
	}
}
