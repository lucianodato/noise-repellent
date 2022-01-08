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

#ifndef NREPEL_H_INCLUDED
#define NREPEL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

// clang-format off
#ifndef NREPEL_EXPORT
#  ifdef _WIN32
#     if defined(NREPEL_BUILD_SHARED) /* build dll */
#         define NREPEL_EXPORT __declspec(dllexport)
#     elif !defined(NREPEL_BUILD_STATIC) /* use dll */
#         define NREPEL_EXPORT __declspec(dllimport)
#     else /* static library */
#         define NREPEL_EXPORT
#     endif
#  else
#     if __GNUC__ >= 4
#         define NREPEL_EXPORT __attribute__((visibility("default")))
#     else
#         define NREPEL_EXPORT
#     endif
#  endif
#endif
// clang-format on

#include <stdbool.h>
#include <stdint.h>

// TODO (luciano/todo): Extract library to it's own repository when API is
// stable
// TODO (luciano/todo): Document interface
// TODO (luciano/todo): Test main file and increase coverage
// TODO (luciano/todo): Move Plugin functionality to the plugin. Circular
// buffer, softbypass, noise_profile declaration, parameters, etc.

typedef void *NoiseRepellentHandle;

typedef struct NrepelDenoiseParameters {
  bool enable;               // Plugin
  bool learn_noise;          // Plugin
  bool residual_listen;      // Plugin
  bool adaptive_noise_learn; // Removed by extracted different plugin
  float reduction_amount;
  float release_time;
  float masking_ceiling_limit; // Should be internal
  float whitening_factor;
  float transient_threshold; // Should be Adaptive
  float noise_rescale;
} NrepelDenoiseParameters;

NREPEL_EXPORT NoiseRepellentHandle nrepel_initialize(uint32_t sample_rate);
NREPEL_EXPORT void nrepel_free(NoiseRepellentHandle instance);
NREPEL_EXPORT bool nrepel_process(NoiseRepellentHandle instance,
                                  uint32_t number_of_samples,
                                  const float *input, float *output);
NREPEL_EXPORT uint32_t nrepel_get_latency(NoiseRepellentHandle instance);
NREPEL_EXPORT uint32_t
nrepel_get_noise_profile_size(NoiseRepellentHandle instance);
NREPEL_EXPORT float *nrepel_get_noise_profile(NoiseRepellentHandle instance);
NREPEL_EXPORT bool nrepel_load_noise_profile(NoiseRepellentHandle instance,
                                             const float *restored_profile,
                                             uint32_t profile_size,
                                             uint32_t profile_blocks);
NREPEL_EXPORT uint32_t
nrepel_get_noise_profile_blocks_averaged(NoiseRepellentHandle instance);
NREPEL_EXPORT bool nrepel_reset_noise_profile(NoiseRepellentHandle instance);
NREPEL_EXPORT bool
nrepel_noise_profile_available(NoiseRepellentHandle instance);
NREPEL_EXPORT bool nrepel_load_parameters(NoiseRepellentHandle instance,
                                          NrepelDenoiseParameters parameters);

// Objective interfaces

// Init and free
// Process adaptive (with parameters)
// Process manual with provided profile (with parameters)
// Estimate noise returns number of blocks averages and out parameters a profile
// (with parameters)

// NREPEL_EXPORT NoiseRepellentHandle nrepel_initialize(int sample_rate, int
// fft_size); NREPEL_EXPORT int nrepel_get_size(NoiseRepellentHandle instance);
// NREPEL_EXPORT void nrepel_free(NoiseRepellentHandle instance);
// NREPEL_EXPORT int nrepel_calculate_noise_profile(NoiseRepellentHandle
// instance,
//                                    const float *input,
//                                    float *calculated_profile);
// NREPEL_EXPORT int nrepel_get_block_size(NoiseRepellentHandle instance);
// NREPEL_EXPORT bool nrepel_denoise_block(NoiseRepellentHandle instance, const
// float *input,
//                           float *output, const float *noise_profile,
//                           float reduction_db, float whitening_percentage,
//                           float release_ms, float noise_gain_db);
// NREPEL_EXPORT bool nrepel_denoise_block_adaptive(NoiseRepellentHandle
// instance,
//                                    const float *input, float *output,
//                                    float reduction_db, float noise_gain_db);

#ifdef __cplusplus
}
#endif
#endif