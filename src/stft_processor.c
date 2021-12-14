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

#include "stft_processor.h"
#include "spectral_utils.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define FFT_SIZE 2048
#define OVERLAP_FACTOR 4
#define INPUT_WINDOW_TYPE 3
#define OUTPUT_WINDOW_TYPE 3

static void stft_processor_pre_and_post_window(STFTProcessor *self);
static void stft_processor_analysis(STFTProcessor *self);
static void stft_processor_synthesis(STFTProcessor *self);

struct STFTProcessor {
  uint32_t fft_size;
  uint32_t half_fft_size;
  fftwf_plan forward;
  fftwf_plan backward;
  uint32_t window_option_input;
  uint32_t window_option_output;
  uint32_t overlap_factor;
  float overlap_scale_factor;
  uint32_t hop;
  uint32_t input_latency;
  uint32_t read_position;
  float *input_window;
  float *output_window;
  float *in_fifo;
  float *out_fifo;
  float *output_accum;
  float *input_fft_buffer;
  float *output_fft_buffer;
};

STFTProcessor *stft_processor_initialize() {
  STFTProcessor *self = (STFTProcessor *)calloc(1, sizeof(STFTProcessor));

  self->fft_size = FFT_SIZE;
  self->half_fft_size = self->fft_size / 2;
  self->window_option_input = INPUT_WINDOW_TYPE;
  self->window_option_output = OUTPUT_WINDOW_TYPE;
  self->overlap_factor = OVERLAP_FACTOR;
  self->hop = self->fft_size / self->overlap_factor;
  self->input_latency = self->fft_size - self->hop;
  self->read_position = self->input_latency;

  self->input_window = (float *)calloc(self->fft_size, sizeof(float));
  self->output_window = (float *)calloc(self->fft_size, sizeof(float));

  self->in_fifo = (float *)calloc(self->fft_size, sizeof(float));
  self->out_fifo = (float *)calloc(self->fft_size, sizeof(float));

  self->output_accum = (float *)calloc((self->fft_size * 2), sizeof(float));

  self->input_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->output_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->forward =
      fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer,
                        self->output_fft_buffer, FFTW_R2HC, FFTW_ESTIMATE);
  self->backward =
      fftwf_plan_r2r_1d(self->fft_size, self->output_fft_buffer,
                        self->input_fft_buffer, FFTW_HC2R, FFTW_ESTIMATE);

  stft_processor_pre_and_post_window(self);

  return self;
}

void stft_processor_free(STFTProcessor *self) {
  free(self->input_fft_buffer);
  free(self->output_fft_buffer);
  fftwf_destroy_plan(self->forward);
  fftwf_destroy_plan(self->backward);
  free(self->input_window);
  free(self->output_window);
  free(self->in_fifo);
  free(self->out_fifo);
  free(self->output_accum);
  free(self);
}

uint32_t get_stft_latency(STFTProcessor *self) { return self->input_latency; }
uint32_t get_fft_size(STFTProcessor *self) { return self->fft_size; }
uint32_t get_overlap_factor(STFTProcessor *self) {
  return self->overlap_factor;
}
uint32_t get_spectral_processing_size(STFTProcessor *self) {
  return self->half_fft_size / 2 + 1;
}

void stft_processor_run(STFTProcessor *self,
                        spectral_processing *spectral_processing,
                        SpectralProcessor *spectral_processor,
                        const uint32_t number_of_samples, const float *input,
                        float *output) {
  for (uint32_t k = 0; k < number_of_samples; k++) {
    self->in_fifo[self->read_position] = input[k];
    output[k] = self->out_fifo[self->read_position - self->input_latency];
    self->read_position++;

    if (self->read_position >= self->fft_size) {
      self->read_position = self->input_latency;

      memcpy(self->input_fft_buffer, self->in_fifo,
             sizeof(float) * self->fft_size);

      stft_processor_analysis(self);

      spectral_processing(spectral_processor, self->output_fft_buffer);

      stft_processor_synthesis(self);
    }
  }
}

static void stft_processor_pre_and_post_window(STFTProcessor *self) {

  switch ((WindowTypes)self->window_option_input) {
  case HANN_WINDOW:
    get_fft_window(self->input_window, self->fft_size,
                   (WindowTypes)self->window_option_input);
    break;
  case HAMMING_WINDOW:
    get_fft_window(self->input_window, self->fft_size,
                   (WindowTypes)self->window_option_input);
    break;
  case BLACKMAN_WINDOW:
    get_fft_window(self->input_window, self->fft_size,
                   (WindowTypes)self->window_option_input);
    break;
  case VORBIS_WINDOW:
    get_fft_window(self->input_window, self->fft_size,
                   (WindowTypes)self->window_option_input);
    break;
  }

  switch ((WindowTypes)self->window_option_output) {
  case HANN_WINDOW:
    get_fft_window(self->output_window, self->fft_size,
                   (WindowTypes)self->window_option_output);
    break;
  case HAMMING_WINDOW:
    get_fft_window(self->output_window, self->fft_size,
                   (WindowTypes)self->window_option_output);
    break;
  case BLACKMAN_WINDOW:
    get_fft_window(self->output_window, self->fft_size,
                   (WindowTypes)self->window_option_output);
    break;
  case VORBIS_WINDOW:
    get_fft_window(self->output_window, self->fft_size,
                   (WindowTypes)self->window_option_output);
    break;
  }

  float sum = 0.f;
  for (uint32_t i = 0; i < self->fft_size; i++) {
    sum += self->input_window[i] * self->output_window[i];
  }

  self->overlap_scale_factor = (sum / (float)(self->fft_size));
}

static void stft_processor_analysis(STFTProcessor *self) {
  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->input_fft_buffer[k] *= self->input_window[k];
  }

  fftwf_execute(self->forward);
}

static void stft_processor_synthesis(STFTProcessor *self) {
  fftwf_execute(self->backward);

  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->input_fft_buffer[k] = self->input_fft_buffer[k] / self->fft_size;
  }

  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->input_fft_buffer[k] =
        (self->output_window[k] * self->input_fft_buffer[k]) /
        (self->overlap_scale_factor * self->overlap_factor);
  }

  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->output_accum[k] += self->input_fft_buffer[k];
  }

  for (uint32_t k = 0; k < self->hop; k++) {
    self->out_fifo[k] = self->output_accum[k];
  }

  memmove(self->output_accum, self->output_accum + self->hop,
          self->fft_size * sizeof(float));

  for (uint32_t k = 0; k < self->input_latency; k++) {
    self->in_fifo[k] = self->in_fifo[k + self->hop];
  }
}