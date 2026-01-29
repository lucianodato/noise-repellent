/*
noise-repellent -- Noise Reduction LV2

Copyright 2022-2026 Luciano Dato <lucianodato@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "../src/noise_profile_state.h"
#include "../src/signal_crossfade.h"
#include "lv2/atom/atom.h"
#include "lv2/core/lv2.h"
#include "lv2/core/lv2_util.h"
#include "lv2/log/logger.h"
#include "lv2/state/state.h"
#include "lv2/urid/urid.h"
#include "specbleach_2d_denoiser.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define NOISEREPELLENT_2D_URI                                                  \
  "https://github.com/lucianodato/noise-repellent#2d"
#define NOISEREPELLENT_2D_STEREO_URI                                           \
  "https://github.com/lucianodato/noise-repellent#2d-stereo"
#define FRAME_SIZE 46

typedef struct URIs {
  LV2_URID atom_Int;
  LV2_URID atom_Float;
  LV2_URID atom_Vector;
  LV2_URID plugin;
  LV2_URID atom_URID;
} URIs;

typedef struct State {
  LV2_URID property_noise_profile_1_mean;
  LV2_URID property_noise_profile_1_median;
  LV2_URID property_noise_profile_1_max;
  LV2_URID property_noise_profile_2_mean;
  LV2_URID property_noise_profile_2_median;
  LV2_URID property_noise_profile_2_max;
  LV2_URID property_noise_profile_size;
  LV2_URID property_averaged_blocks_mean;
  LV2_URID property_averaged_blocks_median;
  LV2_URID property_averaged_blocks_max;
  LV2_URID property_averaged_blocks_2_mean;
  LV2_URID property_averaged_blocks_2_median;
  LV2_URID property_averaged_blocks_2_max;
} State;

static void map_uris(LV2_URID_Map* map, URIs* uris, const char* uri) {
  uris->atom_Int = map->map(map->handle, LV2_ATOM__Int);
  uris->atom_Float = map->map(map->handle, LV2_ATOM__Float);
  uris->atom_Vector = map->map(map->handle, LV2_ATOM__Vector);
  uris->atom_URID = map->map(map->handle, LV2_ATOM__URID);
  uris->plugin = strcmp(uri, NOISEREPELLENT_2D_URI) == 0
                     ? map->map(map->handle, NOISEREPELLENT_2D_URI)
                     : map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI);
}

static void map_state(LV2_URID_Map* map, State* state, const char* uri) {
  if (strcmp(uri, NOISEREPELLENT_2D_STEREO_URI) == 0) {
    state->property_noise_profile_1_mean =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI "#noiseprofilemean");
    state->property_noise_profile_1_median = map->map(
        map->handle, NOISEREPELLENT_2D_STEREO_URI "#noiseprofilemedian");
    state->property_noise_profile_1_max =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI "#noiseprofilemax");
    state->property_noise_profile_2_mean = map->map(
        map->handle, NOISEREPELLENT_2D_STEREO_URI "#noiseprofile2mean");
    state->property_noise_profile_2_median = map->map(
        map->handle, NOISEREPELLENT_2D_STEREO_URI "#noiseprofile2median");
    state->property_noise_profile_2_max =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI "#noiseprofile2max");
    state->property_noise_profile_size =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI "#noiseprofilesize");
    state->property_averaged_blocks_mean =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI
                 "#noiseprofileaveragedblocksmean");
    state->property_averaged_blocks_median =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI
                 "#noiseprofileaveragedblocksmedian");
    state->property_averaged_blocks_max =
        map->map(map->handle,
                 NOISEREPELLENT_2D_STEREO_URI "#noiseprofileaveragedblocksmax");
    state->property_averaged_blocks_2_mean =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI
                 "#noiseprofile2averagedblocksmean");
    state->property_averaged_blocks_2_median =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI
                 "#noiseprofile2averagedblocksmedian");
    state->property_averaged_blocks_2_max =
        map->map(map->handle, NOISEREPELLENT_2D_STEREO_URI
                 "#noiseprofile2averagedblocksmax");

  } else {
    state->property_noise_profile_1_mean =
        map->map(map->handle, NOISEREPELLENT_2D_URI "#noiseprofilemean");
    state->property_noise_profile_1_median =
        map->map(map->handle, NOISEREPELLENT_2D_URI "#noiseprofilemedian");
    state->property_noise_profile_1_max =
        map->map(map->handle, NOISEREPELLENT_2D_URI "#noiseprofilemax");
    state->property_noise_profile_size =
        map->map(map->handle, NOISEREPELLENT_2D_URI "#noiseprofilesize");
    state->property_averaged_blocks_mean = map->map(
        map->handle, NOISEREPELLENT_2D_URI "#noiseprofileaveragedblocksmean");
    state->property_averaged_blocks_median = map->map(
        map->handle, NOISEREPELLENT_2D_URI "#noiseprofileaveragedblocksmedian");
    state->property_averaged_blocks_max = map->map(
        map->handle, NOISEREPELLENT_2D_URI "#noiseprofileaveragedblocksmax");
  }
}

typedef enum PortIndex {
  NOISEREPELLENT_2D_NOISE_LEARN = 0,
  NOISEREPELLENT_2D_MODE = 1,
  NOISEREPELLENT_2D_ADAPTIVE_NOISE = 2,
  NOISEREPELLENT_2D_AMOUNT = 3,
  NOISEREPELLENT_2D_NLM_SMOOTHING = 4,
  NOISEREPELLENT_2D_WHITENING = 5,
  NOISEREPELLENT_2D_RESIDUAL_LISTEN = 6,
  NOISEREPELLENT_2D_BYPASS = 7,
  NOISEREPELLENT_2D_RESET_NOISE_PROFILE = 8,
  NOISEREPELLENT_2D_LATENCY = 9,
  NOISEREPELLENT_2D_INPUT_1 = 10,
  NOISEREPELLENT_2D_OUTPUT_1 = 11,
  NOISEREPELLENT_2D_INPUT_2 = 12,
  NOISEREPELLENT_2D_OUTPUT_2 = 13,
} PortIndex;

typedef struct NoiseRepellent2DPlugin {
  const float* input_1;
  const float* input_2;
  float* input_buf_1;
  float* input_buf_2;
  float* output_1;
  float* output_2;
  float sample_rate;
  float* report_latency;

  LV2_URID_Map* map;
  LV2_Log_Logger log;
  URIs uris;
  State state;
  char* plugin_uri;

  SpectralBleachHandle lib_instance_1;
  SpectralBleachHandle lib_instance_2;
  SpectralBleach2DDenoiserParameters parameters;
  SignalCrossfade* soft_bypass;
  SignalCrossfade* soft_bypass_2;
  NoiseProfileState* noise_profile_state_1;
  NoiseProfileState* noise_profile_state_2;
  float* noise_profile_1;
  float* noise_profile_2;
  uint32_t profile_size;

  float* noise_learn;
  float* noise_reduction_mode;
  float* residual_listen;
  float* reduction_amount;
  float* nlm_smoothing;
  float* whitening;
  float* reset_noise_profile;
  float* bypass;
  float* adaptive_noise;

  bool activated;
  float prev_reset_state;

} NoiseRepellent2DPlugin;

static void cleanup(LV2_Handle instance) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  if (self->noise_profile_state_1) {
    noise_profile_state_free(self->noise_profile_state_1);
    free(self->noise_profile_1);
  }

  if (self->lib_instance_1) {
    specbleach_2d_free(self->lib_instance_1);
  }

  if (self->noise_profile_state_2) {
    noise_profile_state_free(self->noise_profile_state_2);
    free(self->noise_profile_2);
  }

  if (self->lib_instance_2) {
    specbleach_2d_free(self->lib_instance_2);
  }

  if (self->plugin_uri) {
    free(self->plugin_uri);
  }

  if (self->soft_bypass) {
    signal_crossfade_free(self->soft_bypass);
  }
  if (self->soft_bypass_2) {
    signal_crossfade_free(self->soft_bypass_2);
  }

  if (self->input_buf_1) {
    free(self->input_buf_1);
  }

  if (self->input_buf_2) {
    free(self->input_buf_2);
  }

  free(instance);
}

static LV2_Handle instantiate(const LV2_Descriptor* descriptor,
                              const double rate, const char* bundle_path,
                              const LV2_Feature* const* features) {
  NoiseRepellent2DPlugin* self =
      (NoiseRepellent2DPlugin*)calloc(1U, sizeof(NoiseRepellent2DPlugin));
  self->activated = false;
  self->prev_reset_state = 0.F;

  // clang-format off
  const char *missing =
      lv2_features_query(features,
                         LV2_LOG__log, &self->log.log, false,
                         LV2_URID__map, &self->map, true,
                         NULL);
  // clang-format on

  lv2_log_logger_set_map(&self->log, self->map);

  if (missing) {
    lv2_log_error(&self->log, "Missing feature <%s>\n", missing);
    cleanup((LV2_Handle)self);
    return NULL;
  }

  if (strstr(descriptor->URI, NOISEREPELLENT_2D_STEREO_URI)) {
    self->plugin_uri =
        (char*)calloc(strlen(NOISEREPELLENT_2D_STEREO_URI) + 1, sizeof(char));
    strcpy(self->plugin_uri, descriptor->URI);
  } else {
    self->plugin_uri =
        (char*)calloc(strlen(NOISEREPELLENT_2D_URI) + 1, sizeof(char));
    strcpy(self->plugin_uri, descriptor->URI);
  }

  map_uris(self->map, &self->uris, self->plugin_uri);
  map_state(self->map, &self->state, self->plugin_uri);

  self->sample_rate = (float)rate;

  // input_buf_1 is needed for crossfade if the host uses in-place buffers
  self->input_buf_1 = (float*)calloc((size_t)self->sample_rate, sizeof(float));
  if (!self->input_buf_1) {
    lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
    cleanup((LV2_Handle)self);
    return NULL;
  }

  self->lib_instance_1 =
      specbleach_2d_initialize((uint32_t)self->sample_rate, (float)FRAME_SIZE);
  if (!self->lib_instance_1) {
    lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
    cleanup((LV2_Handle)self);
    return NULL;
  }

  self->profile_size =
      specbleach_2d_get_noise_profile_size(self->lib_instance_1);
  self->noise_profile_state_1 =
      noise_profile_state_initialize(self->uris.atom_Float);
  self->noise_profile_1 = (float*)calloc(self->profile_size, sizeof(float));

  if (strstr(self->plugin_uri, NOISEREPELLENT_2D_STEREO_URI)) {
    self->lib_instance_2 = specbleach_2d_initialize((uint32_t)self->sample_rate,
                                                    (float)FRAME_SIZE);
    if (!self->lib_instance_2) {
      lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
      cleanup((LV2_Handle)self);
      return NULL;
    }

    self->noise_profile_state_2 =
        noise_profile_state_initialize(self->uris.atom_Float);
    self->noise_profile_2 = (float*)calloc(self->profile_size, sizeof(float));

    self->input_buf_2 =
        (float*)calloc((size_t)self->sample_rate, sizeof(float));
    if (!self->input_buf_2) {
      lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
      cleanup((LV2_Handle)self);
      return NULL;
    }
  }

  uint32_t latency = specbleach_2d_get_latency(self->lib_instance_1);
  self->soft_bypass =
      signal_crossfade_initialize((uint32_t)self->sample_rate, latency);

  if (!self->soft_bypass) {
    cleanup((LV2_Handle)self);
    return NULL;
  }

  if (self->lib_instance_2) {
    self->soft_bypass_2 =
        signal_crossfade_initialize((uint32_t)self->sample_rate, latency);
    if (!self->soft_bypass_2) {
      cleanup((LV2_Handle)self);
      return NULL;
    }
  }

  return (LV2_Handle)self;
}

static void connect_port(LV2_Handle instance, uint32_t port, void* data) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  switch ((PortIndex)port) {
    case NOISEREPELLENT_2D_NOISE_LEARN:
      self->noise_learn = (float*)data;
      break;
    case NOISEREPELLENT_2D_MODE:
      self->noise_reduction_mode = (float*)data;
      break;
    case NOISEREPELLENT_2D_AMOUNT:
      self->reduction_amount = (float*)data;
      break;
    case NOISEREPELLENT_2D_NLM_SMOOTHING:
      self->nlm_smoothing = (float*)data;
      break;
    case NOISEREPELLENT_2D_WHITENING:
      self->whitening = (float*)data;
      break;
    case NOISEREPELLENT_2D_RESIDUAL_LISTEN:
      self->residual_listen = (float*)data;
      break;
    case NOISEREPELLENT_2D_BYPASS:
      self->bypass = (float*)data;
      break;
    case NOISEREPELLENT_2D_ADAPTIVE_NOISE:
      self->adaptive_noise = (float*)data;
      break;
    case NOISEREPELLENT_2D_RESET_NOISE_PROFILE:
      self->reset_noise_profile = (float*)data;
      break;
    case NOISEREPELLENT_2D_LATENCY:
      self->report_latency = (float*)data;
      break;
    case NOISEREPELLENT_2D_INPUT_1:
      self->input_1 = (const float*)data;
      break;
    case NOISEREPELLENT_2D_OUTPUT_1:
      self->output_1 = (float*)data;
      break;
    default:
      break;
  }
}

static void connect_port_stereo(LV2_Handle instance, uint32_t port,
                                void* data) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  connect_port(instance, port, data);

  switch ((PortIndex)port) {
    case NOISEREPELLENT_2D_INPUT_2:
      self->input_2 = (const float*)data;
      break;
    case NOISEREPELLENT_2D_OUTPUT_2:
      self->output_2 = (float*)data;
      break;
    default:
      break;
  }
}

static void activate(LV2_Handle instance) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  self->activated = true;
  self->prev_reset_state = 0.F;

  if (self->report_latency) {
    *self->report_latency =
        (float)specbleach_2d_get_latency(self->lib_instance_1);
  }

  if (self->soft_bypass) {
    signal_crossfade_reset(self->soft_bypass);
  }
  if (self->soft_bypass_2) {
    signal_crossfade_reset(self->soft_bypass_2);
  }

  if (self->input_buf_1) {
    memset(self->input_buf_1, 0.F, (size_t)self->sample_rate * sizeof(float));
  }

  if (self->input_buf_2) {
    memset(self->input_buf_2, 0.F, (size_t)self->sample_rate * sizeof(float));
  }
}

static void deactivate(LV2_Handle instance) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  self->activated = false;

  if (self->report_latency) {
    *self->report_latency =
        (float)specbleach_2d_get_latency(self->lib_instance_1);
  }
}

static void run(LV2_Handle instance, uint32_t number_of_samples) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  if (!self->input_1 || !self->output_1) {
    return;
  }

  // Handle reset noise profile (edge detection)
  if (self->reset_noise_profile) {
    float current_reset = *self->reset_noise_profile;
    if (current_reset > 0.5F && self->prev_reset_state <= 0.5F) {
      specbleach_2d_reset_noise_profile(self->lib_instance_1);
      if (self->lib_instance_2) {
        specbleach_2d_reset_noise_profile(self->lib_instance_2);
      }
    }
    self->prev_reset_state = current_reset;
  }

  // Copy input for in-place processing
  if (self->input_1 == self->output_1) {
    memcpy(self->input_buf_1, self->input_1,
           sizeof(float) * (number_of_samples <= self->sample_rate
                                ? (size_t)number_of_samples
                                : (size_t)self->sample_rate));
  }

  // clang-format off
  self->parameters = (SpectralBleach2DDenoiserParameters){
      .learn_noise = self->noise_learn ? (int)*self->noise_learn : 0,
      .noise_reduction_mode = self->noise_reduction_mode ? (int)*self->noise_reduction_mode : 1,
      .residual_listen = self->residual_listen ? (bool)*self->residual_listen : false,
      .reduction_amount = self->reduction_amount ? *self->reduction_amount : 10.0f,
      .smoothing_factor = self->nlm_smoothing ? *self->nlm_smoothing : 1.5f,
      .whitening_factor = self->whitening ? *self->whitening : 0.0f,
      .adaptive_noise = self->adaptive_noise ? (int)*self->adaptive_noise : 0,
      .noise_estimation_method = 1, // Always use SPP-MMSE
  };
  // clang-format on

  specbleach_2d_load_parameters(self->lib_instance_1, self->parameters);
  specbleach_2d_process(self->lib_instance_1, number_of_samples, self->input_1,
                        self->output_1);

  bool enable_processing = self->bypass ? ((int)*self->bypass == 0) : true;

  signal_crossfade_run(
      self->soft_bypass, number_of_samples,
      self->input_1 == self->output_1 ? self->input_buf_1 : self->input_1,
      self->output_1, enable_processing);
}

static void run_stereo(LV2_Handle instance, uint32_t number_of_samples) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  run(instance, number_of_samples);

  if (!self->input_2 || !self->output_2 || !self->lib_instance_2) {
    return;
  }

  // Copy input for in-place processing
  if (self->input_2 == self->output_2) {
    memcpy(self->input_buf_2, self->input_2,
           sizeof(float) * (number_of_samples <= self->sample_rate
                                ? (size_t)number_of_samples
                                : (size_t)self->sample_rate));
  }

  specbleach_2d_load_parameters(self->lib_instance_2, self->parameters);
  specbleach_2d_process(self->lib_instance_2, number_of_samples, self->input_2,
                        self->output_2);

  bool enable_processing = self->bypass ? ((int)*self->bypass == 0) : true;

  signal_crossfade_run(
      self->soft_bypass_2, number_of_samples,
      self->input_2 == self->output_2 ? self->input_buf_2 : self->input_2,
      self->output_2, enable_processing);
}

static LV2_State_Status save(LV2_Handle instance,
                             LV2_State_Store_Function store,
                             LV2_State_Handle handle, uint32_t flags,
                             const LV2_Feature* const* features) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  store(handle, self->state.property_noise_profile_size, &self->profile_size,
        sizeof(uint32_t), self->uris.atom_Int,
        LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

  // Save all 3 profiles for channel 1
  for (int mode = 1; mode <= 3; mode++) {
    if (specbleach_2d_noise_profile_available_for_mode(self->lib_instance_1,
                                                       mode)) {
      uint32_t averaged_blocks =
          specbleach_2d_get_noise_profile_blocks_averaged_for_mode(
              self->lib_instance_1, mode);

      LV2_URID blocks_property;
      LV2_URID profile_property;

      if (mode == 1) {
        blocks_property = self->state.property_averaged_blocks_mean;
        profile_property = self->state.property_noise_profile_1_mean;
      } else if (mode == 2) {
        blocks_property = self->state.property_averaged_blocks_median;
        profile_property = self->state.property_noise_profile_1_median;
      } else { // mode == 3
        blocks_property = self->state.property_averaged_blocks_max;
        profile_property = self->state.property_noise_profile_1_max;
      }

      store(handle, blocks_property, &averaged_blocks, sizeof(uint32_t),
            self->uris.atom_Int, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

      memcpy(
          noise_profile_get_elements(self->noise_profile_state_1),
          specbleach_2d_get_noise_profile_for_mode(self->lib_instance_1, mode),
          sizeof(float) * self->profile_size);

      store(handle, profile_property, (void*)self->noise_profile_state_1,
            noise_profile_get_size(), self->uris.atom_Vector,
            LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);
    }
  }

  // Save all 3 profiles for channel 2 (stereo only)
  if (strstr(self->plugin_uri, NOISEREPELLENT_2D_STEREO_URI)) {
    for (int mode = 1; mode <= 3; mode++) {
      if (specbleach_2d_noise_profile_available_for_mode(self->lib_instance_2,
                                                         mode)) {
        uint32_t averaged_blocks =
            specbleach_2d_get_noise_profile_blocks_averaged_for_mode(
                self->lib_instance_2, mode);

        LV2_URID blocks_property;
        LV2_URID profile_property;

        if (mode == 1) {
          blocks_property = self->state.property_averaged_blocks_2_mean;
          profile_property = self->state.property_noise_profile_2_mean;
        } else if (mode == 2) {
          blocks_property = self->state.property_averaged_blocks_2_median;
          profile_property = self->state.property_noise_profile_2_median;
        } else { // mode == 3
          blocks_property = self->state.property_averaged_blocks_2_max;
          profile_property = self->state.property_noise_profile_2_max;
        }

        store(handle, blocks_property, &averaged_blocks, sizeof(uint32_t),
              self->uris.atom_Int, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

        memcpy(noise_profile_get_elements(self->noise_profile_state_2),
               specbleach_2d_get_noise_profile_for_mode(self->lib_instance_2,
                                                        mode),
               sizeof(float) * self->profile_size);

        store(handle, profile_property, (void*)self->noise_profile_state_2,
              noise_profile_get_size(), self->uris.atom_Vector,
              LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);
      }
    }
  }

  return LV2_STATE_SUCCESS;
}

static LV2_State_Status restore(LV2_Handle instance,
                                LV2_State_Retrieve_Function retrieve,
                                LV2_State_Handle handle, uint32_t flags,
                                const LV2_Feature* const* features) {
  NoiseRepellent2DPlugin* self = (NoiseRepellent2DPlugin*)instance;

  size_t size = 0U;
  uint32_t type = 0U;
  uint32_t valflags = 0U;

  const uint32_t* saved_fftsize = (const uint32_t*)retrieve(
      handle, self->state.property_noise_profile_size, &size, &type, &valflags);
  if (saved_fftsize == NULL || type != self->uris.atom_Int) {
    return LV2_STATE_ERR_NO_PROPERTY;
  }
  const uint32_t fftsize = *saved_fftsize;

  // Restore all 3 profiles for channel 1
  for (int mode = 1; mode <= 3; mode++) {
    LV2_URID blocks_property;
    LV2_URID profile_property;

    if (mode == 1) {
      blocks_property = self->state.property_averaged_blocks_mean;
      profile_property = self->state.property_noise_profile_1_mean;
    } else if (mode == 2) {
      blocks_property = self->state.property_averaged_blocks_median;
      profile_property = self->state.property_noise_profile_1_median;
    } else { // mode == 3
      blocks_property = self->state.property_averaged_blocks_max;
      profile_property = self->state.property_noise_profile_1_max;
    }

    const uint32_t* saved_averagedblocks = (const uint32_t*)retrieve(
        handle, blocks_property, &size, &type, &valflags);
    if (saved_averagedblocks == NULL || type != self->uris.atom_Int) {
      continue; // Skip if this profile mode wasn't saved
    }
    const uint32_t averagedblocks = *saved_averagedblocks;

    const void* saved_noise_profile =
        retrieve(handle, profile_property, &size, &type, &valflags);
    if (!saved_noise_profile || size != noise_profile_get_size() ||
        type != self->uris.atom_Vector) {
      continue;
    }

    memcpy(self->noise_profile_1, (float*)LV2_ATOM_BODY(saved_noise_profile),
           sizeof(float) * self->profile_size);

    specbleach_2d_load_noise_profile_for_mode(self->lib_instance_1,
                                              self->noise_profile_1, fftsize,
                                              averagedblocks, mode);
  }

  // Restore all 3 profiles for channel 2 (stereo only)
  if (strstr(self->plugin_uri, NOISEREPELLENT_2D_STEREO_URI)) {
    for (int mode = 1; mode <= 3; mode++) {
      LV2_URID blocks_property;
      LV2_URID profile_property;

      if (mode == 1) {
        blocks_property = self->state.property_averaged_blocks_2_mean;
        profile_property = self->state.property_noise_profile_2_mean;
      } else if (mode == 2) {
        blocks_property = self->state.property_averaged_blocks_2_median;
        profile_property = self->state.property_noise_profile_2_median;
      } else { // mode == 3
        blocks_property = self->state.property_averaged_blocks_2_max;
        profile_property = self->state.property_noise_profile_2_max;
      }

      const uint32_t* saved_averagedblocks = (const uint32_t*)retrieve(
          handle, blocks_property, &size, &type, &valflags);
      if (saved_averagedblocks == NULL || type != self->uris.atom_Int) {
        continue;
      }
      const uint32_t averagedblocks = *saved_averagedblocks;

      const void* saved_noise_profile =
          retrieve(handle, profile_property, &size, &type, &valflags);
      if (!saved_noise_profile || size != noise_profile_get_size() ||
          type != self->uris.atom_Vector) {
        continue;
      }

      memcpy(self->noise_profile_2, (float*)LV2_ATOM_BODY(saved_noise_profile),
             sizeof(float) * self->profile_size);

      specbleach_2d_load_noise_profile_for_mode(self->lib_instance_2,
                                                self->noise_profile_2, fftsize,
                                                averagedblocks, mode);
    }
  }

  return LV2_STATE_SUCCESS;
}

static const void* extension_data(const char* uri) {
  static const LV2_State_Interface state_iface = {save, restore};
  if (!strcmp(uri, LV2_STATE__interface)) {
    return &state_iface;
  }
  return NULL;
}

// clang-format off
static const LV2_Descriptor descriptor_2d = {
    NOISEREPELLENT_2D_URI,
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data
};

static const LV2_Descriptor descriptor_2d_stereo = {
    NOISEREPELLENT_2D_STEREO_URI,
    instantiate,
    connect_port_stereo,
    activate,
    run_stereo,
    deactivate,
    cleanup,
    extension_data
};
// clang-format on

LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor(uint32_t index) {
  switch (index) {
    case 0:
      return &descriptor_2d;
    case 1:
      return &descriptor_2d_stereo;
    default:
      return NULL;
  }
}
