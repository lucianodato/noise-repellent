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

#define PF_SMOOTHING 100.0 //Postfilter smoothing scaling

//General spectral subtraction configuration
#define GAMMA1 2.f
#define GAMMA2 0.5f

//Power substraction
void
power_subtraction(int fft_size_2, float* spectrum, float* noise_thresholds, float* Gk)
{
	int k;

	for (k = 0; k <= fft_size_2 ; k++)
	{
		if (noise_thresholds[k] > FLT_MIN)
		{
			if(spectrum[k] > noise_thresholds[k])
			{
				Gk[k] = (spectrum[k]-noise_thresholds[k]) / spectrum[k];
			}
			else
			{
				Gk[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
			Gk[(2*fft_size_2)-k] = Gk[k];
	}
}

//Gating with envelope smoothing
void
spectral_gating(int fft_size_2, float* spectrum, float* noise_thresholds, float* Gk)
{
	int k;

	for (k = 0; k <= fft_size_2 ; k++)
	{
		if (noise_thresholds[k] > FLT_MIN)
		{
			//Hard knee
			if (spectrum[k] >= noise_thresholds[k])
			{
				//over the threshold
				Gk[k] = 1.f;
			}
			else
			{
				//under the threshold
				Gk[k] = 0.f;
			}
		}
		else
		{
			//Otherwise we keep everything as is
			Gk[k] = 1.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		Gk[(2*fft_size_2)-k] = Gk[k];
	}
}

void
compute_post_filter(int fft_size_2, int fft_size, float* spectrum, float pf_threshold,
										float* postfilter, float* Gk_spectral)
{
	int k;
	float num = 0.f, den = 0.f;
	float indicator;
	float ksi_lambda;
	float n_lambda;

	//Low SNR detector
	for (k = 0; k <= fft_size_2 ; k++)
	{
		num += spectrum[k] * Gk_spectral[k];
		den += spectrum[k];
	}

	indicator = num/den;

	//threshold decision
	if(indicator >= pf_threshold)
	{
		ksi_lambda = 1.f;
	}
	else
	{
		ksi_lambda = indicator;
	}

	//window size
	if(ksi_lambda == 1.f)
	{
		n_lambda = 1.f;
	}
	else
	{
		n_lambda = 2.f*roundf(PF_SMOOTHING*(1.f - ksi_lambda/pf_threshold)) + 1.f;
	}

	//construct the filter window (zero phase)
	for (k = 0; k < fft_size_2 ; k++)
	{
		if(k < n_lambda)
		{
			postfilter[k] = 1.f/n_lambda;
		}
		else
		{
			postfilter[k] = 0.f;
		}
	}

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		postfilter[fft_size-k] = postfilter[k];
	}
}

/*Generalized Spectral Subtraction
 *gamma defines what type of spectral Subtraction is used
 *gamma1=gamma2=1 is magnitude substaction
 *gamma1=2 gamma2=0.5 is power Subtraction
 *gamma1=2 gamma2=1 is wiener filtering
 *alpha is the oversustraction factor
 *beta is the spectral flooring factor
 *reduction_strenght is the other oversustraction designed by the user
 *so there are 2 oversustraction factors
*/

//This version uses an array of alphas and betas
void
denoise_gain_gss(int fft_size_2, float* alpha, float* beta, float* spectrum,
								 float* noise_thresholds, float* Gk)
{
  int k;

  for (k = 0; k <= fft_size_2 ; k++)
	{
    if (spectrum[k] > FLT_MIN)
		{
      if(powf((noise_thresholds[k]/spectrum[k]),GAMMA1) < (1.f/(alpha[k]+beta[k])))
			{
        Gk[k] = MAX(powf(1.f-(alpha[k]*powf((noise_thresholds[k]/spectrum[k]),GAMMA1)),GAMMA2),0.f);
      }
			else
			{
        Gk[k] = MAX(powf(beta[k]*powf((noise_thresholds[k]/spectrum[k]),GAMMA1),GAMMA2),0.f);
      }
    }
		else
		{
      //Otherwise we keep everything as is
      Gk[k] = 1.f;
    }
  }

	//mirrored gain array
	for (k = 1; k < fft_size_2; k++)
	{
		Gk[(2*fft_size_2)-k] = Gk[k];
	}
}
