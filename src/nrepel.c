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

/**
* \file nrepel.c
* \author Luciano Dato
* \brief The main file for host interaction
*/

#include "lv2/lv2plug.in/ns/lv2core/lv2.h"
#include "lv2/lv2plug.in/ns/ext/urid/urid.h"
#include "lv2/lv2plug.in/ns/ext/atom/atom.h"
#include "lv2/lv2plug.in/ns/ext/state/state.h"

#include "spectral_processing.c"

#define NREPEL_URI "https://github.com/lucianodato/noise-repellent"

//STFT default values
#define FFT_SIZE 2048	//Size of the fft transform and frame
#define INPUT_WINDOW 3   //0 HANN 1 HAMMING 2 BLACKMAN 3 VORBIS Input windows for STFT algorithm
#define OUTPUT_WINDOW 3  //0 HANN 1 HAMMING 2 BLACKMAN 3 VORBIS Output windows for STFT algorithm
#define OVERLAP_FACTOR 4 //4 is 75% overlap Values bigger than 4 will rescale correctly (if Vorbis windows is not used)

///---------------------------------------------------------------------

/**
* Enumeration of LV2 ports.
*/
typedef enum {
	NREPEL_AMOUNT = 0,
	NREPEL_NOFFSET = 1,
	NREPEL_RELEASE = 2,
	NREPEL_MASKING = 3,
	NREPEL_T_PROTECT = 4,
	NREPEL_WHITENING = 5,
	NREPEL_N_LEARN = 6,
	NREPEL_N_ADAPTIVE = 7,
	NREPEL_RESET = 8,
	NREPEL_RESIDUAL_LISTEN = 9,
	NREPEL_ENABLE = 10,
	NREPEL_LATENCY = 11,
	NREPEL_INPUT = 12,
	NREPEL_OUTPUT = 13,
} PortIndex;

/**
* Noise Profile state.
*/
typedef struct
{
	uint32_t child_size;
	uint32_t child_type;
	float array[(FFT_SIZE/2) + 1];
} FFTVector;

/**
* Struct for THE noise repellent instance, the host is going to use.
*/
typedef struct
{
	const float *input; //input of samples from host (changing size)
	float *output;		//output of samples to host (changing size)
	float samp_rate;	//Sample rate received from the host

	//Parameters for the algorithm (user input)
	float *amount_of_reduction;		//Amount of noise to reduce in dB
	float *noise_thresholds_offset; //This is to scale the noise profile (over subtraction factor)
	float *release;					//Release time
	float *masking;					//Masking scaling
	float *whitening_factor_pc;		//Whitening amount of the reduction percentage
	float *noise_learn_state;		//Learn Noise state (Manual-Off-Auto)
	float *adaptive_state;			//Autocapture switch
	float *reset_profile;			//Reset Noise switch
	float *residual_listen;			//For noise only listening
	float *transient_protection;	//Multiplier for thresholding onsets with rolling mean
	float *enable;					//For soft bypass (click free bypass)
	float *report_latency;			//Latency necessary

	//Parameters values and arrays for the STFT
	int fft_size;				//FFTW input size
	int fft_size_2;				//FFTW half input size
	int window_option_input;	//Type of input Window for the STFT
	int window_option_output;   //Type of output Window for the STFT
	float overlap_factor;		//oversampling factor for overlap calculations
	float overlap_scale_factor; //Scaling factor for conserving the final amplitude
	int hop;					//Hop size for the STFT
	float *input_window;		//Input Window values
	float *output_window;		//Output Window values

	//Algorithm exta variables
	float tau;						  //time constant for soft bypass
	float wet_dry_target;			  //softbypass target for softbypass
	float wet_dry;					  //softbypass coeff
	float reduction_coeff;			  //Gain to apply to the residual noise
	float release_coeff;			  //Release coefficient for Envelopes
	float amount_of_reduction_linear; //Reduction amount linear value
	float thresholds_offset_linear;   //Threshold offset linear value
	float whitening_factor;			  //Whitening amount of the reduction

	//Buffers for processing and outputting
	int input_latency;
	float *in_fifo;		 //internal input buffer
	float *out_fifo;	 //internal output buffer
	float *output_accum; //FFT output accumulator
	int read_ptr;		 //buffers read pointer

	//FFTW related arrays
	float *input_fft_buffer;
	float *output_fft_buffer;
	fftwf_plan forward;
	fftwf_plan backward;

	//Arrays and variables for getting bins info
	float *fft_p2;		  //power spectrum
	float *fft_magnitude; //magnitude spectrum
	float *fft_phase;	 //phase spectrum

	//noise related
	float *noise_thresholds_p2;		  //captured noise profile power spectrum
	float *noise_thresholds_scaled;   //captured noise profile power spectrum scaled by oversubtraction
	bool noise_thresholds_availables; //indicate whether a noise profile is available or no
	float noise_window_count;		  //Count windows for mean computing

	//smoothing related
	float *smoothed_spectrum;	  //power spectrum to be smoothed
	float *smoothed_spectrum_prev; //previous frame smoothed power spectrum for envelopes

	//Transient preservation related
	float *transient_preserv_prev; //previous frame smoothed power spectrum for envelopes
	float tp_r_mean;
	bool transient_present;
	float tp_window_count;

	//Reduction gains
	float *Gk; //definitive gain

	//Ensemble related
	float *residual_spectrum;
	float *denoised_spectrum;
	float *final_spectrum;

	//whitening related
	float *residual_max_spectrum;
	float max_decay_rate;
	float whitening_window_count;

	//Loizou algorithm
	float *auto_thresholds; //Reference threshold for louizou algorithm
	float *prev_noise_thresholds;
	float *s_pow_spec;
	float *prev_s_pow_spec;
	float *p_min;
	float *prev_p_min;
	float *speech_p_p;
	float *prev_speech_p_p;

	//masking
	float *bark_z;
	float *absolute_thresholds; //absolute threshold of hearing
	float *SSF;
	float *spl_reference_values;
	float *unity_gain_bark_spectrum;
	float *spreaded_unity_gain_bark_spectrum;
	float *alpha_masking;
	float *beta_masking;
	float *input_fft_buffer_at;
	float *output_fft_buffer_at;
	fftwf_plan forward_at;

	//LV2 state URID (Save and restore noise profile)
	LV2_URID_Map *map;
	LV2_URID atom_Vector;
	LV2_URID atom_Int;
	LV2_URID atom_Float;
	LV2_URID prop_fftsize;
	LV2_URID prop_nwindow;
	LV2_URID prop_FFTp2;
} Nrepel;

/**
* Instantiates the plugin.
*/
static LV2_Handle
instantiate(const LV2_Descriptor *descriptor, double rate, const char *bundle_path,
			const LV2_Feature *const *features)
{
	//Actual struct declaration
	Nrepel *self = (Nrepel *)calloc(1, sizeof(Nrepel));

	//Retrieve the URID map callback, and needed URIDs
	for (int i = 0; features[i]; ++i)
	{
		if (!strcmp(features[i]->URI, LV2_URID__map))
		{
			self->map = (LV2_URID_Map *)features[i]->data;
		}
	}
	if (!self->map)
	{
		//bail out: host does not support urid:map
		free(self);
		return NULL;
	}

	//For lv2 state (noise profile saving)
	self->atom_Vector = self->map->map(self->map->handle, LV2_ATOM__Vector);
	self->atom_Int = self->map->map(self->map->handle, LV2_ATOM__Int);
	self->atom_Float = self->map->map(self->map->handle, LV2_ATOM__Float);
	self->prop_fftsize = self->map->map(self->map->handle, NREPEL_URI "#fftsize");
	self->prop_nwindow = self->map->map(self->map->handle, NREPEL_URI "#nwindow");
	self->prop_FFTp2 = self->map->map(self->map->handle, NREPEL_URI "#FFTp2");

	//Sampling related
	self->samp_rate = (float)rate;

	//FFT related
	self->fft_size = FFT_SIZE;
	self->fft_size_2 = self->fft_size / 2;
	self->input_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
	self->output_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
	self->forward = fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer, self->output_fft_buffer, FFTW_R2HC, FFTW_ESTIMATE);
	self->backward = fftwf_plan_r2r_1d(self->fft_size, self->output_fft_buffer, self->input_fft_buffer, FFTW_HC2R, FFTW_ESTIMATE);

	//STFT window related
	self->window_option_input = INPUT_WINDOW;
	self->window_option_output = OUTPUT_WINDOW;
	self->input_window = (float *)calloc(self->fft_size, sizeof(float));
	self->output_window = (float *)calloc(self->fft_size, sizeof(float));

	//buffers for OLA
	self->in_fifo = (float *)calloc(self->fft_size, sizeof(float));
	self->out_fifo = (float *)calloc(self->fft_size, sizeof(float));
	self->output_accum = (float *)calloc(self->fft_size * 2, sizeof(float));
	self->overlap_factor = OVERLAP_FACTOR;
	self->hop = self->fft_size / self->overlap_factor;
	self->input_latency = self->fft_size - self->hop;
	self->read_ptr = self->input_latency; //the initial position because we are that many samples ahead

	//soft bypass
	self->tau = (1.f - expf(-2.f * M_PI * 25.f * 64.f / self->samp_rate));
	self->wet_dry = 0.f;

	//Arrays for getting bins info
	self->fft_p2 = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->fft_magnitude = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->fft_phase = (float *)calloc((self->fft_size_2 + 1), sizeof(float));

	//noise threshold related
	self->noise_thresholds_p2 = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->noise_thresholds_scaled = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->noise_window_count = 0.f;
	self->noise_thresholds_availables = false;

	//noise adaptive estimation related
	self->auto_thresholds = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->prev_noise_thresholds = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->s_pow_spec = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->prev_s_pow_spec = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->p_min = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->prev_p_min = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->speech_p_p = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->prev_speech_p_p = (float *)calloc((self->fft_size_2 + 1), sizeof(float));

	//smoothing related
	self->smoothed_spectrum = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->smoothed_spectrum_prev = (float *)calloc((self->fft_size_2 + 1), sizeof(float));

	//transient preservation
	self->transient_preserv_prev = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->tp_window_count = 0.f;
	self->tp_r_mean = 0.f;
	self->transient_present = false;

	//masking related
	self->bark_z = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->absolute_thresholds = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->unity_gain_bark_spectrum = (float *)calloc(N_BARK_BANDS, sizeof(float));
	self->spreaded_unity_gain_bark_spectrum = (float *)calloc(N_BARK_BANDS, sizeof(float));
	self->spl_reference_values = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->alpha_masking = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->beta_masking = (float *)calloc((self->fft_size_2 + 1), sizeof(float));
	self->SSF = (float *)calloc((N_BARK_BANDS * N_BARK_BANDS), sizeof(float));
	self->input_fft_buffer_at = (float *)calloc(self->fft_size, sizeof(float));
	self->output_fft_buffer_at = (float *)calloc(self->fft_size, sizeof(float));
	self->forward_at = fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer_at, self->output_fft_buffer_at, FFTW_R2HC, FFTW_ESTIMATE);

	//reduction gains related
	self->Gk = (float *)calloc((self->fft_size), sizeof(float));

	//whitening related
	self->residual_max_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->max_decay_rate = expf(-1000.f / (((WHITENING_DECAY_RATE)*self->samp_rate) / self->hop));
	self->whitening_window_count = 0.f;

	//final ensemble related
	self->residual_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->denoised_spectrum = (float *)calloc((self->fft_size), sizeof(float));
	self->final_spectrum = (float *)calloc((self->fft_size), sizeof(float));

	//Window combination initialization (pre processing window post processing window)
	fft_pre_and_post_window(self->input_window, self->output_window,
							self->fft_size, self->window_option_input,
							self->window_option_output, &self->overlap_scale_factor);

	//Set initial gain as unity for the positive part
	initialize_array(self->Gk, 1.f, self->fft_size);

	//Compute adaptive initial thresholds
	compute_auto_thresholds(self->auto_thresholds, self->fft_size, self->fft_size_2,
							self->samp_rate);

	//MASKING initializations
	compute_bark_mapping(self->bark_z, self->fft_size_2, self->samp_rate);
	compute_absolute_thresholds(self->absolute_thresholds, self->fft_size_2,
								self->samp_rate);
	spl_reference(self->spl_reference_values, self->fft_size_2, self->samp_rate,
				  self->input_fft_buffer_at, self->output_fft_buffer_at,
				  &self->forward_at);
	compute_SSF(self->SSF);

	//Initializing unity gain values for offset normalization
	initialize_array(self->unity_gain_bark_spectrum, 1.f, N_BARK_BANDS);
	//Convolve unitary energy bark spectrum with SSF
	convolve_with_SSF(self->SSF, self->unity_gain_bark_spectrum,
					  self->spreaded_unity_gain_bark_spectrum);

	initialize_array(self->alpha_masking, 1.f, self->fft_size_2 + 1);
	initialize_array(self->beta_masking, 0.f, self->fft_size_2 + 1);

	return (LV2_Handle)self;
}

/**
* Used by the host to connect the ports of this plugin.
*/
static void
connect_port(LV2_Handle instance, uint32_t port, void *data)
{
	Nrepel *self = (Nrepel *)instance;

	switch ((PortIndex)port)
	{
	case NREPEL_AMOUNT:
		self->amount_of_reduction = (float *)data;
		break;
	case NREPEL_NOFFSET:
		self->noise_thresholds_offset = (float *)data;
		break;
	case NREPEL_RELEASE:
		self->release = (float *)data;
		break;
	case NREPEL_MASKING:
		self->masking = (float *)data;
		break;
	case NREPEL_WHITENING:
		self->whitening_factor_pc = (float *)data;
		break;
	case NREPEL_N_LEARN:
		self->noise_learn_state = (float *)data;
		break;
	case NREPEL_N_ADAPTIVE:
		self->adaptive_state = (float *)data;
		break;
	case NREPEL_RESIDUAL_LISTEN:
		self->residual_listen = (float *)data;
		break;
	case NREPEL_T_PROTECT:
		self->transient_protection = (float *)data;
		break;
	case NREPEL_RESET:
		self->reset_profile = (float *)data;
		break;
	case NREPEL_ENABLE:
		self->enable = (float *)data;
		break;
	case NREPEL_LATENCY:
		self->report_latency = (float *)data;
		break;
	case NREPEL_INPUT:
		self->input = (const float *)data;
		break;
	case NREPEL_OUTPUT:
		self->output = (float *)data;
		break;
	}
}

/**
* To reset the noise profile and set every value to default one.
*/
static void
reset_noise_profile(Nrepel *self)
{
	initialize_array(self->noise_thresholds_p2, 0.f, self->fft_size_2 + 1);
	initialize_array(self->noise_thresholds_scaled, 0.f, self->fft_size_2 + 1);
	self->noise_window_count = 0.f;
	self->noise_thresholds_availables = false;

	initialize_array(self->Gk, 1.f, self->fft_size);

	initialize_array(self->residual_max_spectrum, 0.f, self->fft_size);
	self->whitening_window_count = 0.f;

	initialize_array(self->prev_noise_thresholds, 0.f, self->fft_size_2 + 1);
	initialize_array(self->s_pow_spec, 0.f, self->fft_size_2 + 1);
	initialize_array(self->prev_s_pow_spec, 0.f, self->fft_size_2 + 1);
	initialize_array(self->p_min, 0.f, self->fft_size_2 + 1);
	initialize_array(self->prev_p_min, 0.f, self->fft_size_2 + 1);
	initialize_array(self->speech_p_p, 0.f, self->fft_size_2 + 1);
	initialize_array(self->prev_speech_p_p, 0.f, self->fft_size_2 + 1);

	initialize_array(self->alpha_masking, 1.f, self->fft_size_2 + 1);
	initialize_array(self->beta_masking, 0.f, self->fft_size_2 + 1);

	self->tp_window_count = 0.f;
	self->tp_r_mean = 0.f;
	self->transient_present = false;
}

/**
* Main process function of the plugin.
*/
static void
run(LV2_Handle instance, uint32_t n_samples)
{
	Nrepel *self = (Nrepel *)instance;

	//handy variables
	int k;
	unsigned int pos;

	//Inform latency at run call
	*(self->report_latency) = (float)self->input_latency;

	//Reset button state (if on)
	if (*(self->reset_profile) == 1.f)
	{
		reset_noise_profile(self);
	}

	//Softbypass targets in case of disabled or enabled
	if (*(self->enable) == 0.f)
	{ //if disabled
		self->wet_dry_target = 0.f;
	}
	else
	{ //if enabled
		self->wet_dry_target = 1.f;
	}
	//Interpolate parameters over time softly to bypass without clicks or pops
	self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry) + FLT_MIN;

	//Parameters values
	/*exponential decay coefficients for envelopes and adaptive noise profiling
		These must take into account the hop size as explained in the following paper
		FFT-BASED DYNAMIC RANGE COMPRESSION*/
	if (*(self->release) != 0.f) //This allows to turn off smoothing with 0 ms in order to use masking only
	{
		self->release_coeff = expf(-1000.f / (((*(self->release)) * self->samp_rate) / self->hop));
	}
	else
	{
		self->release_coeff = 0.f; //This avoids incorrect results when moving sliders rapidly
	}

	self->amount_of_reduction_linear = from_dB(-1.f * *(self->amount_of_reduction));
	self->thresholds_offset_linear = from_dB(*(self->noise_thresholds_offset));
	self->whitening_factor = *(self->whitening_factor_pc) / 100.f;

	//main loop for processing
	for (pos = 0; pos < n_samples; pos++)
	{
		//Store samples int the input buffer
		self->in_fifo[self->read_ptr] = self->input[pos];
		//Output samples in the output buffer (even zeros introduced by latency)
		self->output[pos] = self->out_fifo[self->read_ptr - self->input_latency];
		//Now move the read pointer
		self->read_ptr++;

		//Once the buffer is full we can do stuff
		if (self->read_ptr >= self->fft_size)
		{
			//Reset the input buffer position
			self->read_ptr = self->input_latency;

			//----------STFT Analysis------------

			//Adding and windowing the frame input values in the center (zero-phasing)
			for (k = 0; k < self->fft_size; k++)
			{
				self->input_fft_buffer[k] = self->in_fifo[k] * self->input_window[k];
			}

			//----------FFT Analysis------------

			//Do transform
			fftwf_execute(self->forward);

			//-----------GET INFO FROM BINS--------------

			get_info_from_bins(self->fft_p2, self->fft_magnitude, self->fft_phase,
							   self->fft_size_2, self->fft_size,
							   self->output_fft_buffer);

			/////////////////////SPECTRAL PROCESSING//////////////////////////

			/*This section countains the specific noise reduction processing blocks
				but it could be replaced with any spectral processing (I'm looking at you future tinkerer)
				Parameters for the STFT transform can be changed at the top of this file
			*/

			//If the spectrum is not silence
			if (!is_empty(self->fft_p2, self->fft_size_2))
			{
				//If adaptive noise is selected the noise is adapted in time
				if (*(self->adaptive_state) == 1.f)
				{
					//This has to be revised(issue 8 on github)
					adapt_noise(self->fft_p2, self->fft_size_2, self->noise_thresholds_p2,
								self->auto_thresholds, self->prev_noise_thresholds,
								self->s_pow_spec, self->prev_s_pow_spec, self->p_min,
								self->prev_p_min, self->speech_p_p, self->prev_speech_p_p);

					self->noise_thresholds_availables = true;
				}

				/*If selected estimate noise spectrum is based on selected portion of signal
				 *do not process the signal
				 */
				if (*(self->noise_learn_state) == 1.f)
				{ //MANUAL

					//Increase window count for rolling mean
					self->noise_window_count++;

					get_noise_statistics(self->fft_p2, self->fft_size_2,
										 self->noise_thresholds_p2, self->noise_window_count);

					self->noise_thresholds_availables = true;
				}
				else
				{
					//If there is a noise profile reduce noise
					if (self->noise_thresholds_availables == true)
					{
						//Detector smoothing and oversubtraction
						preprocessing(self->thresholds_offset_linear, self->fft_p2,
									  self->noise_thresholds_p2, self->noise_thresholds_scaled,
									  self->smoothed_spectrum, self->smoothed_spectrum_prev,
									  self->fft_size_2, self->bark_z, self->absolute_thresholds,
									  self->SSF, self->release_coeff,
									  self->spreaded_unity_gain_bark_spectrum,
									  self->spl_reference_values, self->alpha_masking,
									  self->beta_masking, *(self->masking), *(self->adaptive_state),
									  self->amount_of_reduction_linear, self->transient_preserv_prev,
									  &self->tp_window_count, &self->tp_r_mean,
									  &self->transient_present, *(self->transient_protection));

						//Supression rule
						spectral_gain(self->fft_p2, self->noise_thresholds_p2,
									  self->noise_thresholds_scaled, self->smoothed_spectrum,
									  self->fft_size_2, *(self->adaptive_state), self->Gk,
									  *(self->transient_protection), self->transient_present);

						//apply gains
						denoised_calulation(self->fft_size, self->output_fft_buffer,
											self->denoised_spectrum, self->Gk);

						//residual signal
						residual_calulation(self->fft_size, self->output_fft_buffer,
											self->residual_spectrum, self->denoised_spectrum,
											self->whitening_factor, self->residual_max_spectrum,
											&self->whitening_window_count, self->max_decay_rate);

						//Ensemble the final spectrum using residual and denoised
						final_spectrum_ensemble(self->fft_size, self->final_spectrum,
												self->residual_spectrum,
												self->denoised_spectrum,
												self->amount_of_reduction_linear,
												*(self->residual_listen));

						soft_bypass(self->final_spectrum, self->output_fft_buffer,
									self->wet_dry, self->fft_size);
					}
				}
			}

			///////////////////////////////////////////////////////////

			//----------STFT Synthesis------------

			//------------FFT Synthesis-------------

			//Do inverse transform
			fftwf_execute(self->backward);

			//Normalizing value
			for (k = 0; k < self->fft_size; k++)
			{
				self->input_fft_buffer[k] = self->input_fft_buffer[k] / self->fft_size;
			}

			//------------OVERLAPADD-------------

			//Windowing scaling and accumulation
			for (k = 0; k < self->fft_size; k++)
			{
				self->output_accum[k] += (self->output_window[k] * self->input_fft_buffer[k]) / (self->overlap_scale_factor * self->overlap_factor);
			}

			//Output samples up to the hop size
			for (k = 0; k < self->hop; k++)
			{
				self->out_fifo[k] = self->output_accum[k];
			}

			//shift FFT accumulator the hop size
			memmove(self->output_accum, self->output_accum + self->hop,
					self->fft_size * sizeof(float));

			//move input FIFO
			for (k = 0; k < self->input_latency; k++)
			{
				self->in_fifo[k] = self->in_fifo[k + self->hop];
			}
			//-------------------------------
		} //if
	}	 //main loop
}

/**
* Cleanup and freeing memory.
*/
static void
cleanup(LV2_Handle instance)
{
	free(instance);
}

/**
* State saving of the noise profile.
*/
static LV2_State_Status
savestate(LV2_Handle instance, LV2_State_Store_Function store, LV2_State_Handle handle,
		  uint32_t flags, const LV2_Feature *const *features)
{
	Nrepel *self = (Nrepel *)instance;

	FFTVector *vector = (FFTVector *)malloc(sizeof(FFTVector));

	vector->child_type = self->atom_Float;
	vector->child_size = sizeof(float);

	store(handle, self->prop_fftsize, &self->fft_size_2, sizeof(int), self->atom_Int,
		  LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	store(handle, self->prop_nwindow, &self->noise_window_count, sizeof(float),
		  self->atom_Float, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	memcpy(vector->array, self->noise_thresholds_p2, sizeof(vector->array));

	store(handle, self->prop_FFTp2, (void *)vector, sizeof(FFTVector),
		  self->atom_Vector, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

	return LV2_STATE_SUCCESS;
}

/**
* State restoration of the noise profile.
*/
static LV2_State_Status
restorestate(LV2_Handle instance, LV2_State_Retrieve_Function retrieve,
			 LV2_State_Handle handle, uint32_t flags,
			 const LV2_Feature *const *features)
{
	Nrepel *self = (Nrepel *)instance;
	size_t size;
	uint32_t type;
	uint32_t valflags;

	const int32_t *fftsize = retrieve(handle, self->prop_fftsize, &size, &type, &valflags);
	if (!fftsize || type != self->atom_Int || *fftsize != self->fft_size_2)
	{
		return LV2_STATE_ERR_NO_PROPERTY;
	}

	const void *vecFFTp2 = retrieve(handle, self->prop_FFTp2, &size, &type, &valflags);
	if (!vecFFTp2 || size != sizeof(FFTVector) || type != self->atom_Vector)
	{
		return LV2_STATE_ERR_NO_PROPERTY;
	}

	//Deactivate any denoising before loading any noise profile
	self->noise_thresholds_availables = false;

	//Copy to local variables
	memcpy(self->noise_thresholds_p2, (float *)LV2_ATOM_BODY(vecFFTp2), (self->fft_size_2 + 1) * sizeof(float));

	const float *wincount = retrieve(handle, self->prop_nwindow, &size, &type, &valflags);
	if (fftsize && type == self->atom_Float)
	{
		self->noise_window_count = *wincount;
	}

	//Reactivate denoising with restored profile
	self->noise_thresholds_availables = true;

	return LV2_STATE_SUCCESS;
}

/**
* extension for additional interfaces.
*/
static const void *
extension_data(const char *uri)
{
	static const LV2_State_Interface state = {savestate, restorestate};
	if (!strcmp(uri, LV2_STATE__interface))
	{
		return &state;
	}
	return NULL;
}

/**
* Descriptor for linking methods.
*/
static const LV2_Descriptor descriptor =
	{
		NREPEL_URI,
		instantiate,
		connect_port,
		NULL,
		run,
		NULL,
		cleanup,
		extension_data};

/**
* Symbol export using the descriptor above.
*/
LV2_SYMBOL_EXPORT
const LV2_Descriptor *
lv2_descriptor(uint32_t index)
{
	switch (index)
	{
	case 0:
		return &descriptor;
	default:
		return NULL;
	}
}
