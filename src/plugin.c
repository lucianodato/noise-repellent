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

#include "lv2/atom/atom.h"
#include "lv2/core/lv2.h"
#include "lv2/core/lv2_util.h"
#include "lv2/log/logger.h"
#include "lv2/state/state.h"
#include "lv2/urid/urid.h"
#include "noise_repellent.h"
#include <stdlib.h>
#include <string.h>

#define NOISEREPELLENT_URI "https://github.com/lucianodato/noise-repellent"

typedef struct {
  LV2_URID atom_Vector;
  LV2_URID plugin;
  LV2_URID atom_URID;
} URIs;

typedef struct {
  LV2_URID property_saved_noise_profile;
} State;

static inline void map_uris(LV2_URID_Map *map, URIs *uris) {
  uris->plugin = map->map(map->handle, NOISEREPELLENT_URI);
  uris->atom_Vector = map->map(map->handle, LV2_ATOM__Vector);
  uris->atom_URID = map->map(map->handle, LV2_ATOM__URID);
}

static inline void map_state(LV2_URID_Map *map, State *state) {
  state->property_saved_noise_profile =
      map->map(map->handle, NOISEREPELLENT_URI "#savednoiseprofile");
}

typedef enum {
  NOISEREPELLENT_AMOUNT = 0,
  NOISEREPELLENT_NOISE_OFFSET = 1,
  NOISEREPELLENT_RELEASE = 2,
  NOISEREPELLENT_MASKING = 3,
  NOISEREPELLENT_TRANSIENT_PROTECT = 4,
  NOISEREPELLENT_WHITENING = 5,
  NOISEREPELLENT_NOISE_LEARN = 6,
  NOISEREPELLENT_RESIDUAL_LISTEN = 7,
  NOISEREPELLENT_ENABLE = 8,
  NOISEREPELLENT_LATENCY = 9,
  NOISEREPELLENT_INPUT = 10,
  NOISEREPELLENT_OUTPUT = 11,
} PortIndex;

typedef struct {
  const float *input;
  float *output;
  float sample_rate;
  float *report_latency;

  LV2_URID_Map *map;
  LV2_Log_Logger log;
  URIs uris;
  State state;

  NoiseRepellentLib *lib_instance;

  float *enable;
  float *learn_noise;
  float *residual_listen;
  float *reduction_amount;
  float *release_time;
  float *masking_ceiling_limit;
  float *whitening_factor;
  float *transient_threshold;
  float *noise_rescale;

} NoiseRepellent;

static void cleanup(LV2_Handle instance) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  nr_free(self->lib_instance);
  free(instance);
}

static LV2_Handle instantiate(const LV2_Descriptor *descriptor,
                              const double rate, const char *bundle_path,
                              const LV2_Feature *const *features) {
  NoiseRepellent *self = (NoiseRepellent *)calloc(1, sizeof(NoiseRepellent));

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

  map_uris(self->map, &self->uris);
  map_state(self->map, &self->state);

  self->sample_rate = (float)rate;

  self->lib_instance = nr_initialize(self->sample_rate);
  if (!self->lib_instance) {
    lv2_log_error(&self->log, "Error initializing <%s>\n", NOISEREPELLENT_URI);
    cleanup((LV2_Handle)self);
    return NULL;
  }

  return (LV2_Handle)self;
}

static void connect_port(LV2_Handle instance, uint32_t port, void *data) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  switch ((PortIndex)port) {
  case NOISEREPELLENT_AMOUNT:
    self->reduction_amount = (float *)data;
    break;
  case NOISEREPELLENT_NOISE_OFFSET:
    self->noise_rescale = (float *)data;
    break;
  case NOISEREPELLENT_RELEASE:
    self->release_time = (float *)data;
    break;
  case NOISEREPELLENT_MASKING:
    self->masking_ceiling_limit = (float *)data;
    break;
  case NOISEREPELLENT_WHITENING:
    self->whitening_factor = (float *)data;
    break;
  case NOISEREPELLENT_NOISE_LEARN:
    self->learn_noise = (float *)data;
    break;
  case NOISEREPELLENT_RESIDUAL_LISTEN:
    self->residual_listen = (float *)data;
    break;
  case NOISEREPELLENT_TRANSIENT_PROTECT:
    self->transient_threshold = (float *)data;
    break;
  case NOISEREPELLENT_ENABLE:
    self->enable = (float *)data;
    break;
  case NOISEREPELLENT_LATENCY:
    self->report_latency = (float *)data;
    break;
  case NOISEREPELLENT_INPUT:
    self->input = (const float *)data;
    break;
  case NOISEREPELLENT_OUTPUT:
    self->output = (float *)data;
    break;
  }
}

static void run(LV2_Handle instance, uint32_t number_of_samples) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  nr_load_parameters(self->lib_instance, (bool)*self->enable,
                     (bool)*self->learn_noise, *self->masking_ceiling_limit,
                     *self->noise_rescale, *self->reduction_amount,
                     *self->release_time, (bool)*self->residual_listen,
                     *self->transient_threshold, *self->whitening_factor);

  *self->report_latency = (float)nr_get_latency(self->lib_instance);

  nr_process(self->lib_instance, number_of_samples, self->input, self->output);
}

static LV2_State_Status save(LV2_Handle instance,
                             LV2_State_Store_Function store,
                             LV2_State_Handle handle, uint32_t flags,
                             const LV2_Feature *const *features) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  uint32_t noise_profile_size = nr_get_noise_profile_size(self->lib_instance);
  float *noise_profile = nr_get_noise_profile(self->lib_instance);

  store(handle, self->state.property_saved_noise_profile, (void *)noise_profile,
        sizeof(float) * noise_profile_size, self->uris.atom_Vector,
        LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);

  return LV2_STATE_SUCCESS;
}

static LV2_State_Status restore(LV2_Handle instance,
                                LV2_State_Retrieve_Function retrieve,
                                LV2_State_Handle handle, uint32_t flags,
                                const LV2_Feature *const *features) {
  NoiseRepellent *self = (NoiseRepellent *)instance;

  size_t size = 0;
  uint32_t type = 0;
  uint32_t valflags = 0;

  const void *saved_noise_profile =
      retrieve(handle, self->state.property_saved_noise_profile, &size, &type,
               &valflags);
  if (!saved_noise_profile ||
      size != sizeof(float) * nr_get_noise_profile_size(self->lib_instance) ||
      type != self->uris.atom_Vector) {
    return LV2_STATE_ERR_NO_PROPERTY;
  }

  nr_load_noise_profile(self->lib_instance,
                        (float *)LV2_ATOM_BODY(saved_noise_profile));

  return LV2_STATE_SUCCESS;
}

static const void *extension_data(const char *uri) {
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
    NULL,
    run,
    NULL,
    cleanup,
    extension_data
};
// clang-format on

LV2_SYMBOL_EXPORT const LV2_Descriptor *lv2_descriptor(uint32_t index) {
  switch (index) {
  case 0:
    return &descriptor;
  default:
    return NULL;
  }
}
