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

#define MEAN_INTERP 0.5

//------------GAIN AND THRESHOLD CALCULATION---------------

void preprocessing(float noise_thresholds_offset,
							    float* noise_thresholds_scaled,
							    float* smoothed_spectrum,
									float* smoothed_spectrum_prev,
							    int fft_size_2,
							    float* Gk,
									float release_coeff){

	int k;

	//PREPROCESSING

	//------OVERSUBTRACTION------

	//Scale noise thresholds (equals applying an oversubtraction factor in spectral subtraction)
	for (k = 0; k <= fft_size_2; k++) {
		noise_thresholds_scaled[k] *= noise_thresholds_offset;
	}

	//------SMOOTHING DETECTOR------

	/*Time smoothing between current and past Gk (similar effect to ephraim and malah)
		Here is done by applying a release envelope to signal power spectrum
		The best option here is to adaptively smooth 2D spectral components so it will require a biger buffer
		as suggested by Lukin in Suppression of Musical Noise Artifacts in Audio Noise Reduction by Adaptive 2D Filtering
	*/
	apply_envelope(smoothed_spectrum,
								 smoothed_spectrum_prev,
								 fft_size_2,
								 release_coeff);

	memcpy(smoothed_spectrum_prev,smoothed_spectrum,sizeof(float)*(fft_size_2+1));
}

void spectral_gain(float* smoothed_spectrum,
							    float* noise_thresholds_scaled,
							    int fft_size_2,
									float adaptive,
							    float* Gk){

	if(adaptive == 1.f){
		power_subtraction(fft_size_2,
										smoothed_spectrum,
										noise_thresholds_scaled,
										Gk);
	}else{
		spectral_gating(fft_size_2,
										smoothed_spectrum,
										noise_thresholds_scaled,
										Gk);
	}

}

void postprocessing(int fft_size_2,
							      int fft_size,
										float* fft_p2,
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
										float pf_threshold){

  int k;
  float postfilter[fft_size];

	//mirrored gain array to use a power of 2 fft transform
	for (k = 0; k <= fft_size_2; k++) {
		if(k < fft_size_2)
			Gk[fft_size-k] = Gk[k];
	}


	//GAIN SMOOTHING USING A POSTFILTER
	//Compute the filter
	compute_post_filter(fft_size_2,
										fft_size,
										fft_p2,
										pf_threshold,
										postfilter,
										Gk);

	//Convolution using fft transform

	//Copy to fft buffers
	memcpy(input_fft_buffer_ps,postfilter,fft_size*sizeof(float));
	memcpy(input_fft_buffer_g,Gk,fft_size*sizeof(float));

	//FFT Analysis
	fftwf_execute(*forward_ps);
	fftwf_execute(*forward_g);

	//Multiply with the filter computed
	for (k = 0; k <= fft_size_2; k++) {
    output_fft_buffer_g[k] *= output_fft_buffer_ps[k];
    if(k < fft_size_2)
      output_fft_buffer_g[fft_size-k] *= output_fft_buffer_ps[fft_size-k];
  }

	//FFT Synthesis (only gain needs to be synthetised)
	fftwf_execute(*backward_g);

	//Normalizing
	for (k = 0; k < fft_size; k++){
		input_fft_buffer_g[k] = input_fft_buffer_g[k] / fft_size;
	}

	//Copy to orginal arrays
	memcpy(Gk,input_fft_buffer_g,fft_size*sizeof(float));

	///////////////////
}

void denoised_calulation(int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
								      float* denoised_spectrum,
								      float* Gk){

  int k;

  //Apply the computed gain to the signal and store it in denoised array
  for (k = 0; k <= fft_size_2; k++) {
    denoised_spectrum[k] = output_fft_buffer[k] * Gk[k];
    if(k < fft_size_2)
      denoised_spectrum[fft_size-k] = output_fft_buffer[fft_size-k] * Gk[fft_size-k];
  }
}

void residual_calulation(int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
								      float* residual_spectrum,
								      float* denoised_spectrum,
								      float* Gk,
											float whitening_factor){

  int k;

  //Residual signal
  for (k = 0; k <= fft_size_2; k++) {
   residual_spectrum[k] = output_fft_buffer[k] - denoised_spectrum[k];
   if(k < fft_size_2)
    residual_spectrum[fft_size-k] = output_fft_buffer[fft_size-k] - denoised_spectrum[fft_size-k];
  }

	////////////POSTPROCESSING RESIDUAL
	//Whitening (residual spectrum more similar to white noise)
	if(whitening_factor > 0.f) {
		whitening(residual_spectrum,whitening_factor,fft_size_2);
	}
	////////////
}

void final_spectrum_ensemble(int fft_size_2,
								      int fft_size,
								      float* output_fft_buffer,
								      float* residual_spectrum,
								      float* denoised_spectrum,
											float reduction_amount,
								      float wet_dry,
								      float noise_listen){

  int k;

	//OUTPUT RESULTS using smooth bypass and parametric sustraction
	if (noise_listen == 0.f){
	//Mix residual and processed (Parametric way of noise reduction)
		for (k = 0; k <= fft_size_2; k++) {
			output_fft_buffer[k] =  (1.f-wet_dry) * output_fft_buffer[k] + (denoised_spectrum[k] + residual_spectrum[k]*reduction_amount) * wet_dry;
			if(k < fft_size_2)
				output_fft_buffer[fft_size-k] = (1.f-wet_dry) * output_fft_buffer[fft_size-k] + (denoised_spectrum[fft_size-k] + residual_spectrum[fft_size-k]*reduction_amount) * wet_dry;
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
