/*
noise-repellent -- Noise Reduction LV2

Copyright 2021 Luciano Dato <lucianodato@gmail.com>

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

#include "../src/signal_crossfade.h"
#include "../subprojects/libnrepel/include/nrepel.h"
#include "lv2/core/lv2.h"
#include <stdlib.h>

#define NOISEREPELLENT_ADAPTIVE_URI                                            \
  "https://github.com/lucianodato/noise-repellent#adaptive"
#define NOISEREPELLENT_ADAPTIVE_STEREO_URI                                     \
  "https://github.com/lucianodato/noise-repellent#adaptive-stereo"

// TODO (luciano/todo): Use state mapping and unmapping instead of ladspa float
// arguments

typedef enum PortIndex {
  NOISEREPELLENT_AMOUNT = 0,
  NOISEREPELLENT_NOISE_OFFSET = 1,
  NOISEREPELLENT_RESIDUAL_LISTEN = 2,
  NOISEREPELLENT_ENABLE = 3,
  NOISEREPELLENT_LATENCY = 4,
  NOISEREPELLENT_INPUT_1 = 5,
  NOISEREPELLENT_OUTPUT_1 = 6,
  NOISEREPELLENT_INPUT_2 = 7,
  NOISEREPELLENT_OUTPUT_2 = 8,
} PortIndex;

typedef struct NoiseRepellentAdaptivePlugin {
  const float *input_1;
  const float *input_2;
  float *output_1;
  float *output_2;
  float sample_rate;
  float *report_latency;

  NoiseRepellentHandle lib_instance;
  NrepelDenoiseParameters parameters;
  SignalCrossfade *soft_bypass;

  float *enable;
  float *residual_listen;
  float *reduction_amount;
  float *noise_rescale;

} NoiseRepellentAdaptivePlugin;

static void cleanup_adaptive(LV2_Handle instance) {
  NoiseRepellentAdaptivePlugin *self = (NoiseRepellentAdaptivePlugin *)instance;

  if (self->lib_instance) {
    nrepel_free(self->lib_instance);
  }
  if (self->soft_bypass) {
    signal_crossfade_free(self->soft_bypass);
  }

  free(instance);
}

static LV2_Handle instantiate_adaptive(const LV2_Descriptor *descriptor,
                                       const double rate,
                                       const char *bundle_path,
                                       const LV2_Feature *const *features) {
  NoiseRepellentAdaptivePlugin *self = (NoiseRepellentAdaptivePlugin *)calloc(
      1U, sizeof(NoiseRepellentAdaptivePlugin));

  self->sample_rate = (float)rate;

  self->lib_instance = nrepel_initialize((uint32_t)self->sample_rate);
  if (!self->lib_instance) {
    cleanup_adaptive((LV2_Handle)self);
    return NULL;
  }

  self->soft_bypass = signal_crossfade_initialize(self->sample_rate);

  if (!self->soft_bypass) {
    nrepel_free(self);
    return NULL;
  }

  return (LV2_Handle)self;
}

static void connect_port_adaptive(LV2_Handle instance, uint32_t port,
                                  void *data) {
  NoiseRepellentAdaptivePlugin *self = (NoiseRepellentAdaptivePlugin *)instance;

  switch ((PortIndex)port) {
  case NOISEREPELLENT_AMOUNT:
    self->reduction_amount = (float *)data;
    break;
  case NOISEREPELLENT_NOISE_OFFSET:
    self->noise_rescale = (float *)data;
    break;
  case NOISEREPELLENT_RESIDUAL_LISTEN:
    self->residual_listen = (float *)data;
    break;
  case NOISEREPELLENT_ENABLE:
    self->enable = (float *)data;
    break;
  case NOISEREPELLENT_LATENCY:
    self->report_latency = (float *)data;
    break;
  case NOISEREPELLENT_INPUT_1:
    self->input_1 = (const float *)data;
    break;
  case NOISEREPELLENT_OUTPUT_1:
    self->output_1 = (float *)data;
    break;
  default:
    break;
  }
}

static void connect_port_adaptive_stereo(LV2_Handle instance, uint32_t port,
                                         void *data) {
  NoiseRepellentAdaptivePlugin *self = (NoiseRepellentAdaptivePlugin *)instance;

  connect_port_adaptive(instance, port, data);

  switch ((PortIndex)port) {
  case NOISEREPELLENT_INPUT_2:
    self->input_2 = (const float *)data;
    break;
  case NOISEREPELLENT_OUTPUT_2:
    self->output_2 = (float *)data;
    break;
  default:
    break;
  }
}

static void activate_adaptive(LV2_Handle instance) {
  NoiseRepellentAdaptivePlugin *self = (NoiseRepellentAdaptivePlugin *)instance;

  *self->report_latency = (float)nrepel_get_latency(self->lib_instance);
}

static void run_adaptive(LV2_Handle instance, uint32_t number_of_samples) {
  NoiseRepellentAdaptivePlugin *self = (NoiseRepellentAdaptivePlugin *)instance;

  // clang-format off
  self->parameters = (NrepelDenoiseParameters){
      .residual_listen = (bool)*self->residual_listen,
      .reduction_amount = *self->reduction_amount,
      .noise_rescale = *self->noise_rescale
  };
  // clang-format on

  nrepel_load_parameters(self->lib_instance, self->parameters);

  nrepel_process(self->lib_instance, number_of_samples, self->input_1,
                 self->output_1);

  signal_crossfade_run(self->soft_bypass, number_of_samples, self->input_1,
                       self->output_1, (bool)*self->enable);
}

static void run_adaptive_stereo(LV2_Handle instance,
                                uint32_t number_of_samples) {
  NoiseRepellentAdaptivePlugin *self = (NoiseRepellentAdaptivePlugin *)instance;

  run_adaptive(instance, number_of_samples); // Call left side first

  nrepel_process(self->lib_instance, number_of_samples, self->input_2,
                 self->output_2);

  signal_crossfade_run(self->soft_bypass, number_of_samples, self->input_2,
                       self->output_2, (bool)*self->enable);
}

static const void *extension_data(const char *uri) { return NULL; }

// clang-format off
static const LV2_Descriptor descriptor_adaptive = {
    NOISEREPELLENT_ADAPTIVE_URI,
    instantiate_adaptive,
    connect_port_adaptive,
    activate_adaptive,
    run_adaptive,
    NULL,
    cleanup_adaptive,
    extension_data
};

static const LV2_Descriptor descriptor_adaptive_stereo = {
    NOISEREPELLENT_ADAPTIVE_STEREO_URI,
    instantiate_adaptive,
    connect_port_adaptive_stereo,
    activate_adaptive,
    run_adaptive_stereo,
    NULL,
    cleanup_adaptive,
    extension_data
};
// clang-format on

LV2_SYMBOL_EXPORT const LV2_Descriptor *lv2_descriptor(uint32_t index) {
  switch (index) {
  case 0:
    return &descriptor_adaptive;
  case 1:
    return &descriptor_adaptive_stereo;
  default:
    return NULL;
  }
}
