/*
libspecbleach - A spectral processing library

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

#ifndef SPECBLEACH_H_INCLUDED
#define SPECBLEACH_H_INCLUDED

#include <stdbool.h>
#include <stdint.h>

// TODO (luciano/todo): Extract library to it's own repository when API is
// stable. Manage visibility with meson
// TODO (luciano/todo): Document interface
// TODO (luciano/todo): Test main file and increase coverage
// TODO (luciano/todo): Add post filter for frequency smoothing
// TODO (luciano/todo): Move Plugin functionality to the plugin.
// noise_profile instance, parameters, etc.

typedef void *SpectralBleachHandle;

typedef struct SpectralBleachParameters {
  bool learn_noise;     // Plugin
  bool residual_listen; // Plugin
  float reduction_amount;
  float release_time;
  float masking_ceiling_limit; // Should be internal
  float whitening_factor;
  float transient_threshold; // Should be Adaptive or fixed
  float noise_rescale;
} SpectralBleachParameters;

SpectralBleachHandle specbleach_initialize(uint32_t sample_rate);
void specbleach_free(SpectralBleachHandle instance);
bool specbleach_process(SpectralBleachHandle instance,
                        uint32_t number_of_samples, const float *input,
                        float *output);
uint32_t specbleach_get_latency(SpectralBleachHandle instance);
bool specbleach_load_parameters(SpectralBleachHandle instance,
                                SpectralBleachParameters parameters);
uint32_t specbleach_get_noise_profile_size(SpectralBleachHandle instance);
float *specbleach_get_noise_profile(SpectralBleachHandle instance);
bool specbleach_load_noise_profile(SpectralBleachHandle instance,
                                   const float *restored_profile,
                                   uint32_t profile_size,
                                   uint32_t profile_blocks);
uint32_t
specbleach_get_noise_profile_blocks_averaged(SpectralBleachHandle instance);
bool specbleach_reset_noise_profile(SpectralBleachHandle instance);
bool specbleach_noise_profile_available(SpectralBleachHandle instance);

SpectralBleachHandle specbleach_adaptive_initialize(uint32_t sample_rate);
void specbleach_adaptive_free(SpectralBleachHandle instance);
bool specbleach_adaptive_process(SpectralBleachHandle instance,
                                 uint32_t number_of_samples, const float *input,
                                 float *output);
uint32_t specbleach_adaptive_get_latency(SpectralBleachHandle instance);
bool specbleach_adaptive_load_parameters(SpectralBleachHandle instance,
                                         SpectralBleachParameters parameters);

#endif