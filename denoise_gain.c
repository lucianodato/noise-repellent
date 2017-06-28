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
	float gain, Fk;

	for (k = 0; k <= fft_size_2 ; k++) {
		if (noise_thresholds[k] > FLT_MIN){
			if(spectrum[k] > 0.f){
				gain = MAX(spectrum[k]-noise_thresholds[k], 0.f) / spectrum[k];
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

//Gating with envelope smoothing
void gating(int fft_size_2,
	    float fs,
	    float* spectrum,
	    float* noise_thresholds,
	    float* Gk_gate,
	    float* Gk_prev_gate) {

	int k;

	float attack = expf(-logf(9.f)/(fs*0.01));//1ms
	float release = expf(-logf(9.f)/(fs*0.05));//50ms

	for (k = 0; k <= fft_size_2 ; k++) {
		if (noise_thresholds[k] > FLT_MIN){
			//Envelopes application
			if (spectrum[k] >= noise_thresholds[k]){
				Gk_gate[k] = 1.f; // only avoid applying reduction if over the threshold
			}else{
				Gk_gate[k] = 0.f;
			}

			if (Gk_gate[k] > Gk_prev_gate[k])
				Gk_gate[k] = (1.f-attack)*Gk_prev_gate[k] + attack*Gk_gate[k];
			else
				Gk_gate[k] = (1.f-release)*Gk_prev_gate[k] + release*Gk_gate[k];

			//update previous gain
			Gk_prev_gate[k] = Gk_gate[k];
		} else {
			//Otherwise we keep everything as is
			Gk_gate[k] = 1.f;
		}
	}
}

//Both power sustraction and gating at the same time
void hybrid_reduction(int fft_size_2,
									    float fs,
									    float* spectrum,
									    float* noise_thresholds,
									    float* Gk,
									    float* Gk_prev,
											float gsmoothing) {

	int k;

	float gain, Fk, trigger, alpha;
	float attack = expf(-logf(9.f)/(fs*0.01));//1ms
	float release = expf(-logf(9.f)/(fs*0.05));//50ms

	for (k = 0; k <= fft_size_2 ; k++) {
		if (noise_thresholds[k] > FLT_MIN){
			//Spectral sustraction for the frequency
			if(spectrum[k] > 0.f){
				alpha = sqrtf(spectrum[k]/noise_thresholds[k]);
				gain = MAX(spectrum[k] - alpha*noise_thresholds[k], 0.f) / spectrum[k];
			}else{
				gain = 0.f;
			}

			//Avoid invalid gain numbers
			Fk = (1.f-gain);
			if(Fk < 0.f) Fk = 0.f;
			if(Fk > 1.f) Fk = 1.f;
			gain =  1.f - Fk;

			//Gate triggering
			if(spectrum[k] >= noise_thresholds[k]){
				trigger = 1.f;
			} else {
				trigger = 0.f;
			}

			//Hybrid gain
			Gk[k] = gsmoothing*trigger + (1.f-gsmoothing)*gain;

			//Applying envelopes
			if (Gk[k] > Gk_prev[k]){
				Gk[k] = (1.f-attack)*Gk_prev[k] + attack*Gk[k];
			}else{
				Gk[k] = (1.f-release)*Gk_prev[k] + release*Gk[k];
			}

			//update previous gain
			Gk_prev[k] = Gk[k];

		} else {
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}
}
