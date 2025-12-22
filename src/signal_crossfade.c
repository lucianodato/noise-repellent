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

#include "signal_crossfade.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535F
#endif

#define RELEASE_TIME_MS 30.F

struct SignalCrossfade {
  float tau;
  float wet_dry_target;
  float wet_dry;
  float* delay_buffer;
  uint32_t latency;
  uint32_t rw_ptr;
  uint32_t samples_processed;
};

SignalCrossfade* signal_crossfade_initialize(const uint32_t sample_rate,
                                             const uint32_t latency) {
  SignalCrossfade* self = (SignalCrossfade*)calloc(1U, sizeof(SignalCrossfade));

  if (!self) {
    return NULL;
  }

  self->tau =
      (1.F - expf(-128.F * M_PI * RELEASE_TIME_MS / (float)sample_rate));
  self->wet_dry = 0.F;
  self->wet_dry_target = 0.F;
  self->latency = latency;
  self->samples_processed = 0;

  if (latency > 0) {
    self->delay_buffer = (float*)calloc(latency, sizeof(float));
    if (!self->delay_buffer) {
      free(self);
      return NULL;
    }
  }

  return self;
}

static void signal_crossfade_update_wetdry_target(SignalCrossfade* self,
                                                  const bool enable) {
  if (enable) {
    self->wet_dry_target = 1.F;
  } else {
    self->wet_dry_target = 0.F;
  }

  self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry);
}

void signal_crossfade_reset(SignalCrossfade* self) {
  if (!self) {
    return;
  }

  self->wet_dry = 0.F;
  self->wet_dry_target = 0.F;
  self->rw_ptr = 0;
  self->samples_processed = 0;

  // Clear the delay buffer
  if (self->delay_buffer && self->latency > 0) {
    memset(self->delay_buffer, 0, self->latency * sizeof(float));
  }
}

void signal_crossfade_free(SignalCrossfade* self) {
  if (self) {
    if (self->delay_buffer) {
      free(self->delay_buffer);
    }
    free(self);
  }
}


bool signal_crossfade_run(SignalCrossfade* self,
                          const uint32_t number_of_samples, const float* input,
                          float* output, const bool enable) {
  if (!input || !output || number_of_samples <= 0U) {
    return false;
  }

  signal_crossfade_update_wetdry_target(self, enable);

  // Always apply latency compensation to ensure dry and wet signals are aligned
  for (uint32_t k = 0U; k < number_of_samples; k++) {
    float delayed_dry;
    float current_wet = output[k];  // The processed output (already in output buffer)

    // During initial latency period, pass through to avoid startup artifacts
    if (self->samples_processed < self->latency) {
      // Still fill the delay buffer for future use
      if (self->delay_buffer) {
        self->delay_buffer[self->rw_ptr] = input[k];
        self->rw_ptr++;
        if (self->rw_ptr >= self->latency) {
          self->rw_ptr = 0;
        }
      }
      delayed_dry = input[k];  // Use current input as "delayed" during startup
    } else {
      // After latency period, use properly delayed input
      if (self->delay_buffer) {
        // Get the delayed input from the buffer
        delayed_dry = self->delay_buffer[self->rw_ptr];
        // Store current input for future delay
        self->delay_buffer[self->rw_ptr] = input[k];
        self->rw_ptr++;
        if (self->rw_ptr >= self->latency) {
          self->rw_ptr = 0;
        }
      } else {
        delayed_dry = input[k];  // No latency
      }
    }

    // Crossfade between aligned delayed dry and delayed wet signals
    output[k] = (1.F - self->wet_dry) * delayed_dry + self->wet_dry * current_wet;

    self->samples_processed++;
  }

  return true;
}
