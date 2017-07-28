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

#define G_THRESH 0.5 //Fixed threshold to differenciate low and high SNR zones (this should be estimated automatically)

<<<<<<< HEAD
//Power subtraction
=======
#define G_THRESH 0.5    //wideband gate threshold for the detector

//Power substraction
>>>>>>> 1ca2b97... Implemented PostFilter. Not yet working but mostly there. Artifact control and smoothing are kind of a joke soundwise. Release should be delayed to include this functionality
void power_subtraction(int fft_size_2,
		       float* spectrum,
		       float* noise_thresholds,
		       float* Gk) {

	int k;

	for (k = 0; k <= fft_size_2 ; k++) {
		if (noise_thresholds[k] > FLT_MIN){
			if(spectrum[k] > noise_thresholds[k]){
				Gk[k] = (spectrum[k]-noise_thresholds[k]) / spectrum[k];
			} else {
				Gk[k] = 0.f;
			}
		} else {
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}
}

//Gating with envelope smoothing
void spectral_gating(int fft_size_2,
	    float* spectrum,
	    float* noise_thresholds,
	    float* Gk) {

	int k;

	for (k = 0; k <= fft_size_2 ; k++) {
		if (noise_thresholds[k] > FLT_MIN){
			//Hard knee
			if (spectrum[k] >= noise_thresholds[k]){
				//over the threshold
				Gk[k] = 1.f;
			}else{
				//under the threshold
				Gk[k] = 0.f;
			}
		} else {
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}
}

//This probably could be better (Postfilter) TODO
float wideband_gating(int fft_size_2,
	    float* spectrum,
	    float* noise_thresholds,
			float* Gk) {

	int k;
	float num = 0.f, den = 0.f;
	float indicator;
	float Gk_gate;

	//Using low level SNR DETECTOR
	for (k = 0; k <= fft_size_2 ; k++) {
		num += spectrum[k]*Gk[k];
		den += spectrum[k];
	}

	indicator = num/den;

	//Hard knee decision (this has to be somewhat adaptive) TODO
	if (indicator >= G_THRESH){
		//over the threshold
		Gk_gate = 1.f;
	}else{
		//under the threshold
		Gk_gate = 0.f;
	}

	return Gk_gate;
}

void compute_post_filter(int fft_size_2,
									int fft_size,
									float* spectrum,
									float snr_threshold,
									float amount_smoothing,
									float* postfilter,
									float* Gk_spectral) {

	int k;
	float num = 0.f, den = 0.f;
	float indicator;
	float ksi_lambda;
	float n_lambda;

	//Low SNR detector
	for (k = 0; k <= fft_size_2 ; k++) {
		num += spectrum[k] * Gk_spectral[k];
		den += spectrum[k];
		// if(k < fft_size_2){
		// 	num += spectrum[fft_size - k] * Gk_spectral[fft_size - k];
		// 	den += spectrum[fft_size - k];
		// }
	}

	indicator = num/den;

	//threshold decision
	if(indicator >= snr_threshold){
		ksi_lambda = 1.f;
	}else{
		ksi_lambda = indicator;
	}

	//soft decision
	if(ksi_lambda == 1.f){
		n_lambda = 1.f;
	}else{
		n_lambda = 2.f*roundf(amount_smoothing*(1.f - ksi_lambda/snr_threshold)) + 1.f;
	}

	//construct the filter window
	for (k = 0; k <= fft_size_2 ; k++) {
		postfilter[k] = 1.f;
		// if(k < n_lambda){
		// 	postfilter[k] = 1.f/n_lambda;
		// }else{
		// 	postfilter[k] = 0.f;
		// }
		if(k < fft_size_2){
				postfilter[fft_size - k] = postfilter[k];//Mirrored value
		}
	}
}
