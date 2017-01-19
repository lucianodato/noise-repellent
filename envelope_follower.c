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

// //attack and release coefficients
// nrepel->attack_coeff = powf( 0.01, 1.0 / ( *(nrepel->attack) * nrepel->samp_rate) * 0.001 );
// nrepel->release_coeff = powf( 0.01, 1.0 / ( *(nrepel->release) * nrepel->samp_rate) * 0.001 );

// //Envelope Follower
// for (k = 0; k <= nrepel->fft_size_2; k++) {
// 	if(nrepel->fft_magnitude[k] > nrepel->envelope[k]){
// 		nrepel->envelope[k] =  nrepel->attack_coeff*(nrepel->envelope[k] - nrepel->fft_magnitude[k]) + nrepel->fft_magnitude[k];
// 	}else{
// 		nrepel->envelope[k] = nrepel->release_coeff*(nrepel->envelope[k] - nrepel->fft_magnitude[k]) + nrepel->fft_magnitude[k];
// 	}
// 	//Attenuate signal below the threshold following the envelope
// 	if(nrepel->envelope[k] <= nrepel->noise_thresholds[k]){
// 		nrepel->Gk[k] -= nrepel->envelope[k];
// 	}
// }
