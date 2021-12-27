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

#include "signal_crossfade.h"
#include "../shared/configurations.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct SignalCrossfade {
  float tau;
  float wet_dry_target;
  float wet_dry;
};

SignalCrossfade *signal_crossfade_initialize(const uint32_t sample_rate) {
  SignalCrossfade *self =
      (SignalCrossfade *)calloc(1U, sizeof(SignalCrossfade));

  self->tau = (1.F - expf(-2.F * M_PI * 25.F * 64.F / (float)sample_rate));
  self->wet_dry = 0.F;

  return self;
}

void signal_crossfade_free(SignalCrossfade *self) { free(self); }

static void signal_crossfade__update_wetdry_target(SignalCrossfade *self,
                                                   const bool enable) {
  if (enable) {
    self->wet_dry_target = 1.F;
  } else {
    self->wet_dry_target = 0.F;
  }

  self->wet_dry += self->tau * (self->wet_dry_target - self->wet_dry) + FLT_MIN;
}

bool signal_crossfade_run(SignalCrossfade *self,
                          const uint32_t number_of_samples, const float *input,
                          float *output, const bool enable) {
  if (!input || !output || number_of_samples <= 0) {
    return false;
  }

  signal_crossfade__update_wetdry_target(self, enable);

  for (uint32_t k = 0U; k < number_of_samples; k++) {
    output[k] = (1.F - self->wet_dry) * input[k] + output[k] * self->wet_dry;
  }

  return true;
}