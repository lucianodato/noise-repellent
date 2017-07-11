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

//Non linear Power Sustraction
void nonlinear_power_sustraction(float snr_influence,
				 int fft_size_2,
				 float* spectrum,
				 float* noise_thresholds,
				 float* Gk) {
	int k;
	float gain, Fk, alpha;

	for (k = 0; k <= fft_size_2 ; k++) {
		if (noise_thresholds[k] > FLT_MIN){
			if(spectrum[k] > 0.f){
				if(snr_influence > 0.f){
					alpha = snr_influence + sqrtf(spectrum[k]/noise_thresholds[k]);
				}else{
					alpha = 1.f;//Non linear spectral sustraction off
				}
				gain = MAX(spectrum[k]-alpha*noise_thresholds[k], 0.f) / spectrum[k];
			} else {
				gain = 0.f;
			}

			//Avoid invalid gain numbers
			Fk = (1.f-gain);

			if(Fk < 0.f) Fk = 0.f;
			if(Fk > 1.f) Fk = 1.f;

			Gk[k] =  1.f - Fk;

		} else {
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}
}

//Power Sustraction
void power_sustraction(int fft_size_2,
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
void gating(int fft_size_2,
	    float attack_coeff,
			float release_coeff,
	    float* spectrum,
	    float* noise_thresholds,
	    float* Gk,
	    float* Gk_prev) {

	int k;
	float current_value,knee_width;

	for (k = 0; k <= fft_size_2 ; k++) {
		if (noise_thresholds[k] > FLT_MIN){
			//Envelopes application
			if (spectrum[k] >= noise_thresholds[k]){
				current_value = 1.f; // only avoid applying reduction if over the threshold
			}else{
				current_value = 0.f;
			}

			//Applying envelopes
			if (current_value < Gk_prev[k])
				//Is starting to reduce
				Gk[k] = attack_coeff*Gk_prev[k] + (1.f-attack_coeff)*current_value;
			else
				//Is finishing reduction
				Gk[k] = release_coeff*Gk_prev[k] + (1.f-release_coeff)*current_value;
		} else {
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
		//update previous gain
		Gk_prev[k] = Gk[k];
	}
}
