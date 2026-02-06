/*
noise-repellent -- Noise Reduction LV2

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

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
#include "specbleach_denoiser.h"
#include <stdlib.h>
#include <string.h>

#define NOISEREPELLENT_URI "https://github.com/lucianodato/noise-repellent#new"
#define NOISEREPELLENT_STEREO_URI                                              \
  "https://github.com/lucianodato/noise-repellent-stereo#new"
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
  uris->plugin = strcmp(uri, NOISEREPELLENT_URI) == 0
                     ? map->map(map->handle, NOISEREPELLENT_URI)
                     : map->map(map->handle, NOISEREPELLENT_STEREO_URI);
  uris->atom_Int = map->map(map->handle, LV2_ATOM__Int);
  uris->atom_Float = map->map(map->handle, LV2_ATOM__Float);
  uris->atom_Vector = map->map(map->handle, LV2_ATOM__Vector);
  uris->atom_URID = map->map(map->handle, LV2_ATOM__URID);
}

static void map_state(LV2_URID_Map* map, State* state, const char* uri) {
  if (strcmp(uri, NOISEREPELLENT_STEREO_URI) == 0) {
    state->property_noise_profile_1_mean =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI "#noiseprofilemean");
    state->property_noise_profile_1_median =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI "#noiseprofilemedian");
    state->property_noise_profile_1_max =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI "#noiseprofilemax");
    state->property_noise_profile_2_mean =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI "#noiseprofile2mean");
    state->property_noise_profile_2_median =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI "#noiseprofile2median");
    state->property_noise_profile_2_max =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI "#noiseprofile2max");
    state->property_noise_profile_size =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI "#noiseprofilesize");
    state->property_averaged_blocks_mean =
        map->map(map->handle,
                 NOISEREPELLENT_STEREO_URI "#noiseprofileaveragedblocksmean");
    state->property_averaged_blocks_median =
        map->map(map->handle,
                 NOISEREPELLENT_STEREO_URI "#noiseprofileaveragedblocksmedian");
    state->property_averaged_blocks_max =
        map->map(map->handle,
                 NOISEREPELLENT_STEREO_URI "#noiseprofileaveragedblocksmax");
    state->property_averaged_blocks_2_mean =
        map->map(map->handle,
                 NOISEREPELLENT_STEREO_URI "#noiseprofile2averagedblocksmean");
    state->property_averaged_blocks_2_median =
        map->map(map->handle, NOISEREPELLENT_STEREO_URI
                 "#noiseprofile2averagedblocksmedian");
    state->property_averaged_blocks_2_max =
        map->map(map->handle,
                 NOISEREPELLENT_STEREO_URI "#noiseprofile2averagedblocksmax");

  } else {
    state->property_noise_profile_1_mean =
        map->map(map->handle, NOISEREPELLENT_URI "#noiseprofilemean");
    state->property_noise_profile_1_median =
        map->map(map->handle, NOISEREPELLENT_URI "#noiseprofilemedian");
    state->property_noise_profile_1_max =
        map->map(map->handle, NOISEREPELLENT_URI "#noiseprofilemax");
    state->property_noise_profile_size =
        map->map(map->handle, NOISEREPELLENT_URI "#noiseprofilesize");
    state->property_averaged_blocks_mean = map->map(
        map->handle, NOISEREPELLENT_URI "#noiseprofileaveragedblocksmean");
    state->property_averaged_blocks_median = map->map(
        map->handle, NOISEREPELLENT_URI "#noiseprofileaveragedblocksmedian");
    state->property_averaged_blocks_max = map->map(
        map->handle, NOISEREPELLENT_URI "#noiseprofileaveragedblocksmax");
  }
}

typedef enum PortIndex {
  NOISEREPELLENT_NOISE_LEARN = 0,
  NOISEREPELLENT_MODE = 1,
  NOISEREPELLENT_ADAPTIVE_NOISE = 2,
  NOISEREPELLENT_ADAPTIVE_METHOD = 3,
  NOISEREPELLENT_AMOUNT = 4,
  NOISEREPELLENT_NOISE_REDUCTION_TYPE = 5,
  NOISEREPELLENT_NOISE_OFFSET = 6,
  NOISEREPELLENT_SMOOTHING = 7,
  NOISEREPELLENT_WHITENING = 8,
  NOISEREPELLENT_RESIDUAL_LISTEN = 9,
  NOISEREPELLENT_RESET_NOISE_PROFILE = 10,
  NOISEREPELLENT_BYPASS = 11,
  NOISEREPELLENT_LATENCY = 12,
  NOISEREPELLENT_INPUT_1 = 13,
  NOISEREPELLENT_OUTPUT_1 = 14,
  NOISEREPELLENT_INPUT_2 = 15,
  NOISEREPELLENT_OUTPUT_2 = 16,
} PortIndex;

typedef struct NoiseRepellentPlugin {
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

  SignalCrossfade* soft_bypass;
  SignalCrossfade* soft_bypass_2;
  SpectralBleachHandle lib_instance_1;
  SpectralBleachHandle lib_instance_2;
  SpectralBleachDenoiserParameters parameters;
  NoiseProfileState* noise_profile_state_1;
  NoiseProfileState* noise_profile_state_2;
  float* noise_profile_1;
  float* noise_profile_2;
  uint32_t profile_size;

  float* learn_noise;
  float* mode;
  float* noise_scaling_type;
  float* residual_listen;
  float* reduction_amount;
  float* smoothing_factor;
  float* whitening_factor;
  float* noise_rescale;
  float* reset_noise_profile;
  float* bypass;
  float* adaptive_noise;
  float* adaptive_method;

  bool activated;

} NoiseRepellentPlugin;

static void cleanup(LV2_Handle instance) {
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  if (self->noise_profile_state_1) {
    noise_profile_state_free(self->noise_profile_state_1);
    free(self->noise_profile_1);
  }

  if (self->input_buf_1) {
    free(self->input_buf_1);
  }

  if (self->input_buf_2) {
    free(self->input_buf_2);
  }

  if (self->lib_instance_1) {
    specbleach_free(self->lib_instance_1);
  }

  if (self->noise_profile_state_2) {
    noise_profile_state_free(self->noise_profile_state_2);
    free(self->noise_profile_2);
  }

  if (self->lib_instance_2) {
    specbleach_free(self->lib_instance_2);
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

  free(instance);
}

static LV2_Handle instantiate(const LV2_Descriptor* descriptor,
                              const double rate, const char* bundle_path,
                              const LV2_Feature* const* features) {
  NoiseRepellentPlugin* self =
      (NoiseRepellentPlugin*)calloc(1U, sizeof(NoiseRepellentPlugin));
  self->activated = false;

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

  if (strstr(descriptor->URI, NOISEREPELLENT_STEREO_URI)) {
    self->plugin_uri =
        (char*)calloc(strlen(NOISEREPELLENT_STEREO_URI) + 1U, sizeof(char));
    strcpy(self->plugin_uri, descriptor->URI);
  } else {
    self->plugin_uri =
        (char*)calloc(strlen(NOISEREPELLENT_URI) + 1U, sizeof(char));
    strcpy(self->plugin_uri, descriptor->URI);
  }

  map_uris(self->map, &self->uris, self->plugin_uri);
  map_state(self->map, &self->state, self->plugin_uri);

  self->sample_rate = (float)rate;

  // input_buf_1 is needed for crossfade if the host uses in-place buffers
  // (e.g. Ardour)
  self->input_buf_1 = (float*)calloc((size_t)self->sample_rate, sizeof(float));
  if (!self->input_buf_1) {
    lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
    cleanup((LV2_Handle)self);
    return NULL;
  }

  self->lib_instance_1 =
      specbleach_initialize((uint32_t)self->sample_rate, (float)FRAME_SIZE);
  if (!self->lib_instance_1) {
    lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
    cleanup((LV2_Handle)self);
    return NULL;
  }

  self->profile_size = specbleach_get_noise_profile_size(self->lib_instance_1);
  lv2_log_note(&self->log, "Saved Noise Repellent Profile Size <%u>\n",
               (unsigned int)self->profile_size);
  self->noise_profile_state_1 =
      noise_profile_state_initialize(self->uris.atom_Float);

  self->noise_profile_1 = (float*)calloc(self->profile_size, sizeof(float));

  if (strstr(self->plugin_uri, NOISEREPELLENT_STEREO_URI)) {
    self->lib_instance_2 =
        specbleach_initialize((uint32_t)self->sample_rate, (float)FRAME_SIZE);

    if (!self->lib_instance_2) {
      lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
      cleanup((LV2_Handle)self);
      return NULL;
    }

    self->noise_profile_state_2 =
        noise_profile_state_initialize(self->uris.atom_Float);

    self->noise_profile_2 = (float*)calloc(self->profile_size, sizeof(float));

    // input_buf_2 is needed for crossfade if the host uses in-place buffers
    // (e.g. Ardour)
    self->input_buf_2 =
        (float*)calloc((size_t)self->sample_rate, sizeof(float));
    if (!self->input_buf_2) {
      lv2_log_error(&self->log, "Error initializing <%s>\n", self->plugin_uri);
      cleanup((LV2_Handle)self);
      return NULL;
    }
  }

  uint32_t latency = specbleach_get_latency(self->lib_instance_1);
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
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  switch ((PortIndex)port) {
    case NOISEREPELLENT_NOISE_LEARN:
      self->learn_noise = (float*)data;
      break;
    case NOISEREPELLENT_AMOUNT:
      self->reduction_amount = (float*)data;
      break;
    case NOISEREPELLENT_NOISE_REDUCTION_TYPE:
      self->noise_scaling_type = (float*)data;
      break;
    case NOISEREPELLENT_NOISE_OFFSET:
      self->noise_rescale = (float*)data;
      break;
    case NOISEREPELLENT_SMOOTHING:
      self->smoothing_factor = (float*)data;
      break;
    case NOISEREPELLENT_WHITENING:
      self->whitening_factor = (float*)data;
      break;
    case NOISEREPELLENT_MODE:
      self->mode = (float*)data;
      break;
    case NOISEREPELLENT_RESET_NOISE_PROFILE:
      self->reset_noise_profile = (float*)data;
      break;
    case NOISEREPELLENT_BYPASS:
      self->bypass = (float*)data;
      break;
    case NOISEREPELLENT_ADAPTIVE_NOISE:
      self->adaptive_noise = (float*)data;
      break;
    case NOISEREPELLENT_ADAPTIVE_METHOD:
      self->adaptive_method = (float*)data;
      break;
    case NOISEREPELLENT_RESIDUAL_LISTEN:
      self->residual_listen = (float*)data;
      break;
    case NOISEREPELLENT_LATENCY:
      self->report_latency = (float*)data;
      break;
    case NOISEREPELLENT_INPUT_1:
      self->input_1 = (const float*)data;
      break;
    case NOISEREPELLENT_OUTPUT_1:
      self->output_1 = (float*)data;
      break;
    default:
      break;
  }
}

static void connect_port_stereo(LV2_Handle instance, uint32_t port,
                                void* data) {
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  connect_port(instance, port, data);

  switch ((PortIndex)port) {
    case NOISEREPELLENT_INPUT_2:
      self->input_2 = (const float*)data;
      break;
    case NOISEREPELLENT_OUTPUT_2:
      self->output_2 = (float*)data;
      break;
    default:
      break;
  }
}

static void activate(LV2_Handle instance) {
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  self->activated = true;
  // Always report the correct processing latency when activated
  if (self->report_latency) {
    *self->report_latency = (float)specbleach_get_latency(self->lib_instance_1);
  }

  // Reset soft bypass to dry state for clean reactivation
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
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  self->activated = false;
  // Keep reporting the same latency for consistency - soft bypass handles
  // alignment
  if (self->report_latency) {
    *self->report_latency = (float)specbleach_get_latency(self->lib_instance_1);
  }
}

static void run(LV2_Handle instance, uint32_t number_of_samples) {
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  // Check for valid audio ports
  if (!self->input_1 || !self->output_1) {
    return;
  }

  // If the host gave us in-place buffers, we need to keep a copy of the input
  if (self->input_1 == self->output_1) {
    memcpy(self->input_buf_1, self->input_1,
           sizeof(float) * (number_of_samples <= self->sample_rate
                                ? (size_t)number_of_samples
                                : (size_t)self->sample_rate));
  }

  // clang-format off
  self->parameters = (SpectralBleachDenoiserParameters){
      .learn_noise = self->learn_noise ? (int)*self->learn_noise : 0,
      .noise_reduction_mode = self->mode ? (int)*self->mode : 1,
      .residual_listen = self->residual_listen ? (bool)*self->residual_listen : false,
      .noise_scaling_type = self->noise_scaling_type ? (int)*self->noise_scaling_type : 2,
      .reduction_amount = self->reduction_amount ? *self->reduction_amount : 10.0f,
      .noise_rescale = self->noise_rescale ? *self->noise_rescale : 2.0f,
      .smoothing_factor = self->smoothing_factor ? *self->smoothing_factor : 0.0f,
      .whitening_factor = self->whitening_factor ? *self->whitening_factor : 0.0f,
      .adaptive_noise = self->adaptive_noise ? (int)*self->adaptive_noise : 0,
      .noise_estimation_method = self->adaptive_method ? (int)*self->adaptive_method : 0,
  };
  // clang-format on

  specbleach_load_parameters(self->lib_instance_1, self->parameters);

  if (self->reset_noise_profile && (bool)*self->reset_noise_profile) {
    specbleach_reset_noise_profile(self->lib_instance_1);
  }

  // Always process for both learning and denoising
  specbleach_process(self->lib_instance_1, number_of_samples, self->input_1,
                     self->output_1);

  // Apply soft bypass (wet/dry crossfade)
  bool enable_processing = self->bypass ? ((int)*self->bypass == 0) : true;

  signal_crossfade_run(
      self->soft_bypass, number_of_samples,
      self->input_1 == self->output_1 ? self->input_buf_1 : self->input_1,
      self->output_1, enable_processing);
}

static void run_stereo(LV2_Handle instance, uint32_t number_of_samples) {
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  run(instance, number_of_samples);

  // If the host gave us in-place buffers, we need to keep a copy of the input
  if (self->input_2 == self->output_2) {
    memcpy(self->input_buf_2, self->input_2,
           sizeof(float) * (number_of_samples <= self->sample_rate
                                ? (size_t)number_of_samples
                                : (size_t)self->sample_rate));
  }

  specbleach_load_parameters(self->lib_instance_2, self->parameters);

  if (self->reset_noise_profile && (bool)*self->reset_noise_profile) {
    specbleach_reset_noise_profile(self->lib_instance_2);
  }

  // Always process for both learning and denoising
  specbleach_process(self->lib_instance_2, number_of_samples, self->input_2,
                     self->output_2);

  // Apply soft bypass based on bypass parameter and learning mode
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
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

  store(handle, self->state.property_noise_profile_size, &self->profile_size,
        sizeof(uint32_t), self->uris.atom_Int,
        LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

  // Save all 3 profiles for channel 1
  for (int mode = 1; mode <= 3; mode++) {
    if (specbleach_noise_profile_available_for_mode(self->lib_instance_1,
                                                    mode)) {
      uint32_t averaged_blocks =
          specbleach_get_noise_profile_blocks_averaged_for_mode(
              self->lib_instance_1, mode);

      LV2_URID blocks_property;
      LV2_URID profile_property;

      if (mode == 1) {
        blocks_property = self->state.property_averaged_blocks_mean;
        profile_property = self->state.property_noise_profile_1_mean;
      } else if (mode == 2) {
        blocks_property = self->state.property_averaged_blocks_median;
        profile_property = self->state.property_noise_profile_1_median;
      } else {
        blocks_property = self->state.property_averaged_blocks_max;
        profile_property = self->state.property_noise_profile_1_max;
      }

      store(handle, blocks_property, &averaged_blocks, sizeof(uint32_t),
            self->uris.atom_Int, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

      memcpy(noise_profile_get_elements(self->noise_profile_state_1),
             specbleach_get_noise_profile_for_mode(self->lib_instance_1, mode),
             sizeof(float) * self->profile_size);

      store(handle, profile_property, (void*)self->noise_profile_state_1,
            noise_profile_get_size(), self->uris.atom_Vector,
            LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);
    }
  }

  // Save all 3 profiles for channel 2 (stereo only)
  if (strstr(self->plugin_uri, NOISEREPELLENT_STEREO_URI)) {
    for (int mode = 1; mode <= 3; mode++) {
      if (specbleach_noise_profile_available_for_mode(self->lib_instance_2,
                                                      mode)) {
        uint32_t averaged_blocks =
            specbleach_get_noise_profile_blocks_averaged_for_mode(
                self->lib_instance_2, mode);

        LV2_URID blocks_property;
        LV2_URID profile_property;

        if (mode == 1) {
          blocks_property = self->state.property_averaged_blocks_2_mean;
          profile_property = self->state.property_noise_profile_2_mean;
        } else if (mode == 2) {
          blocks_property = self->state.property_averaged_blocks_2_median;
          profile_property = self->state.property_noise_profile_2_median;
        } else {
          blocks_property = self->state.property_averaged_blocks_2_max;
          profile_property = self->state.property_noise_profile_2_max;
        }

        store(handle, blocks_property, &averaged_blocks, sizeof(uint32_t),
              self->uris.atom_Int, LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

        memcpy(
            noise_profile_get_elements(self->noise_profile_state_2),
            specbleach_get_noise_profile_for_mode(self->lib_instance_2, mode),
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
  NoiseRepellentPlugin* self = (NoiseRepellentPlugin*)instance;

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
    } else {
      blocks_property = self->state.property_averaged_blocks_max;
      profile_property = self->state.property_noise_profile_1_max;
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

    memcpy(self->noise_profile_1, (float*)LV2_ATOM_BODY(saved_noise_profile),
           sizeof(float) * self->profile_size);

    specbleach_load_noise_profile_for_mode(self->lib_instance_1,
                                           self->noise_profile_1, fftsize,
                                           averagedblocks, mode);
  }

  // Restore all 3 profiles for channel 2 (stereo only)
  if (strstr(self->plugin_uri, NOISEREPELLENT_STEREO_URI)) {
    for (int mode = 1; mode <= 3; mode++) {
      LV2_URID blocks_property;
      LV2_URID profile_property;

      if (mode == 1) {
        blocks_property = self->state.property_averaged_blocks_2_mean;
        profile_property = self->state.property_noise_profile_2_mean;
      } else if (mode == 2) {
        blocks_property = self->state.property_averaged_blocks_2_median;
        profile_property = self->state.property_noise_profile_2_median;
      } else {
        blocks_property = self->state.property_averaged_blocks_2_max;
        profile_property = self->state.property_noise_profile_2_max;
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
        continue; // Skip if this profile mode wasn't saved
      }

      memcpy(self->noise_profile_2, (float*)LV2_ATOM_BODY(saved_noise_profile),
             sizeof(float) * self->profile_size);

      specbleach_load_noise_profile_for_mode(self->lib_instance_2,
                                             self->noise_profile_2, fftsize,
                                             averagedblocks, mode);
    }
  }

  return LV2_STATE_SUCCESS;
}

static const void* extension_data(const char* uri) {
  static const LV2_State_Interface state = {save, restore};
  if (strcmp(uri, LV2_STATE__interface) == 0) {
    return &state;
  }
  return NULL;
}

// clang-format off
static const LV2_Descriptor descriptor = {
    NOISEREPELLENT_URI,
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data
};

static const LV2_Descriptor descriptor_stereo = {
    NOISEREPELLENT_STEREO_URI,
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
      return &descriptor;
    case 1:
      return &descriptor_stereo;
    default:
      return NULL;
  }
}
