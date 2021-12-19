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
#include "../shared/spectral_utils.h"
#include "stft_windows.h"
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define FFT_SIZE 2048
#define OVERLAP_FACTOR 4

static void get_overlap_scale_factor(StftProcessor *self);
static void stft_analysis(StftProcessor *self);
static void stft_synthesis(StftProcessor *self);
static void stft_write_buffer(StftProcessor *self);
void static stft_transform_and_process(StftProcessor *self,
                                       spectral_processing *spectral_processing,
                                       SPECTAL_PROCESSOR spectral_processor);

typedef struct {
  float overlap_scale_factor;
  uint32_t input_latency;
  uint32_t read_position;
  uint32_t hop;
  uint32_t overlap_factor;
  float *in_fifo;
  float *out_fifo;
  float *output_accumulator;
} SamplesBuffer;

typedef struct {
  float *power_spectrum;
  float *phase_spectrum;
  float *magnitude_spectrum;
} ProcessingSpectrums;

struct StftProcessor {
  uint32_t fft_size;
  uint32_t half_fft_size;
  fftwf_plan forward;
  fftwf_plan backward;

  float *input_fft_buffer;
  float *output_fft_buffer;

  SamplesBuffer stft_buffer;
  StftWindows *stft_windows;

  ProcessingSpectrums processing_spectrums;
};

StftProcessor *stft_processor_initialize() {
  StftProcessor *self = (StftProcessor *)calloc(1, sizeof(StftProcessor));

  self->fft_size = FFT_SIZE;
  self->half_fft_size = self->fft_size / 2;

  self->stft_buffer.overlap_factor = OVERLAP_FACTOR;
  self->stft_buffer.hop = self->fft_size / self->stft_buffer.overlap_factor;
  self->stft_buffer.input_latency = self->fft_size - self->stft_buffer.hop;
  self->stft_buffer.read_position = self->stft_buffer.input_latency;
  self->stft_buffer.in_fifo = (float *)calloc(self->fft_size, sizeof(float));
  self->stft_buffer.out_fifo = (float *)calloc(self->fft_size, sizeof(float));
  self->stft_buffer.output_accumulator =
      (float *)calloc((self->fft_size * 2), sizeof(float));

  self->input_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->output_fft_buffer = (float *)calloc(self->fft_size, sizeof(float));
  self->forward =
      fftwf_plan_r2r_1d(self->fft_size, self->input_fft_buffer,
                        self->output_fft_buffer, FFTW_R2HC, FFTW_ESTIMATE);
  self->backward =
      fftwf_plan_r2r_1d(self->fft_size, self->output_fft_buffer,
                        self->input_fft_buffer, FFTW_HC2R, FFTW_ESTIMATE);

  self->processing_spectrums.power_spectrum =
      (float *)calloc((self->half_fft_size + 1), sizeof(float));

  self->stft_windows = stft_window_initialize(self->fft_size);

  get_overlap_scale_factor(self);

  return self;
}

void stft_processor_free(StftProcessor *self) {
  free(self->input_fft_buffer);
  free(self->output_fft_buffer);
  fftwf_destroy_plan(self->forward);
  fftwf_destroy_plan(self->backward);

  free(self->processing_spectrums.power_spectrum);

  free(self->stft_buffer.in_fifo);
  free(self->stft_buffer.out_fifo);
  free(self->stft_buffer.output_accumulator);

  stft_window_free(self->stft_windows);

  free(self);
}

uint32_t get_stft_latency(StftProcessor *self) {
  return self->stft_buffer.input_latency;
}
uint32_t get_fft_size(StftProcessor *self) { return self->fft_size; }
uint32_t get_overlap_factor(StftProcessor *self) {
  return self->stft_buffer.overlap_factor;
}
uint32_t get_spectral_processing_size(StftProcessor *self) {
  return self->half_fft_size / 2 + 1;
}

static void stft_read_buffer(StftProcessor *self,
                             const uint32_t number_of_samples,
                             const float *input, float *output) {
  for (uint32_t k = 0; k < number_of_samples; k++) {
    self->stft_buffer.in_fifo[self->stft_buffer.read_position] = input[k];
    output[k] = self->stft_buffer.out_fifo[self->stft_buffer.read_position -
                                           self->stft_buffer.input_latency];
    self->stft_buffer.read_position++;

    if (self->stft_buffer.read_position >= self->fft_size) {
      break;
    }
  }
}

void stft_processor_run(StftProcessor *self,
                        spectral_processing *spectral_processing,
                        SPECTAL_PROCESSOR spectral_processor,
                        const uint32_t number_of_samples, const float *input,
                        float *output) {
  for (uint32_t k = 0; k < number_of_samples; k++) {
    self->stft_buffer.in_fifo[self->stft_buffer.read_position] = input[k];
    output[k] = self->stft_buffer.out_fifo[self->stft_buffer.read_position -
                                           self->stft_buffer.input_latency];
    self->stft_buffer.read_position++;

    if (self->stft_buffer.read_position >= self->fft_size) {
      self->stft_buffer.read_position = self->stft_buffer.input_latency;

      memcpy(self->input_fft_buffer, self->stft_buffer.in_fifo,
             sizeof(float) * self->fft_size);

      stft_transform_and_process(self, spectral_processing, spectral_processor);

      stft_write_buffer(self);
    }
  }
}

void static stft_transform_and_process(StftProcessor *self,
                                       spectral_processing *spectral_processing,
                                       SPECTAL_PROCESSOR spectral_processor) {
  stft_analysis(self);

  get_fft_power_spectrum(self->output_fft_buffer, self->fft_size,
                         self->processing_spectrums.power_spectrum,
                         self->half_fft_size);

  spectral_processing(spectral_processor,
                      self->processing_spectrums.power_spectrum);

  memcpy(self->output_fft_buffer, self->processing_spectrums.power_spectrum,
         sizeof(float) * self->half_fft_size + 1);

  stft_synthesis(self);
}

static void get_overlap_scale_factor(StftProcessor *self) {
  self->stft_buffer.overlap_scale_factor =
      (input_output_window_sum(self->stft_windows) / (float)(self->fft_size)) *
      self->stft_buffer.overlap_factor;
}

static void stft_analysis(StftProcessor *self) {
  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->input_fft_buffer[k] *= get_input_window(self->stft_windows)[k];
  }

  fftwf_execute(self->forward);
}

static void stft_synthesis(StftProcessor *self) {
  fftwf_execute(self->backward);

  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->input_fft_buffer[k] =
        (get_output_window(self->stft_windows)[k] * self->input_fft_buffer[k]) /
        (self->stft_buffer.overlap_scale_factor * self->fft_size);
  }
}

static void stft_write_buffer(StftProcessor *self) {
  for (uint32_t k = 0; k < self->fft_size; k++) {
    self->stft_buffer.output_accumulator[k] += self->input_fft_buffer[k];
  }

  for (uint32_t k = 0; k < self->stft_buffer.hop; k++) {
    self->stft_buffer.out_fifo[k] = self->stft_buffer.output_accumulator[k];
  }

  memmove(self->stft_buffer.output_accumulator,
          self->stft_buffer.output_accumulator + self->stft_buffer.hop,
          self->fft_size * sizeof(float));

  for (uint32_t k = 0; k < self->stft_buffer.input_latency; k++) {
    self->stft_buffer.in_fifo[k] =
        self->stft_buffer.in_fifo[k + self->stft_buffer.hop];
  }
}