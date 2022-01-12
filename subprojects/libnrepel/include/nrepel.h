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

#include <stdbool.h>
#include <stdint.h>

// TODO (luciano/todo): Extract library to it's own repository when API is
// stable. Manage visibility with meson
// TODO (luciano/todo): Document interface
// TODO (luciano/todo): Test main file and increase coverage
// TODO (luciano/todo): Move Plugin functionality to the plugin. Softbypass,
// noise_profile declaration, parameters, etc.

typedef void *NoiseRepellentHandle;

typedef struct NrepelDenoiseParameters {
  bool enable;          // Plugin
  bool learn_noise;     // Plugin
  bool residual_listen; // Plugin
  float reduction_amount;
  float release_time;
  float masking_ceiling_limit; // Should be internal
  float whitening_factor;
  float transient_threshold; // Should be Adaptive or fixed
  float noise_rescale;
} NrepelDenoiseParameters;

NoiseRepellentHandle nrepel_initialize(uint32_t sample_rate);
void nrepel_free(NoiseRepellentHandle instance);
bool nrepel_process(NoiseRepellentHandle instance, uint32_t number_of_samples,
                    const float *input, float *output);
bool nrepel_process_adaptive(NoiseRepellentHandle instance,
                             uint32_t number_of_samples, const float *input,
                             float *output);
uint32_t nrepel_get_latency(NoiseRepellentHandle instance);
uint32_t nrepel_get_noise_profile_size(NoiseRepellentHandle instance);
float *nrepel_get_noise_profile(NoiseRepellentHandle instance);
bool nrepel_load_noise_profile(NoiseRepellentHandle instance,
                               const float *restored_profile,
                               uint32_t profile_size, uint32_t profile_blocks);
uint32_t
nrepel_get_noise_profile_blocks_averaged(NoiseRepellentHandle instance);
bool nrepel_reset_noise_profile(NoiseRepellentHandle instance);
bool nrepel_noise_profile_available(NoiseRepellentHandle instance);
bool nrepel_load_parameters(NoiseRepellentHandle instance,
                            NrepelDenoiseParameters parameters);

// Objective interfaces

// Init and free
// Process adaptive (with parameters)
// Process manual with provided profile (with parameters)
// Estimate noise returns number of blocks averages and out parameters a profile
// (with parameters)

//  NoiseRepellentHandle nrepel_initialize(int sample_rate, int
// fft_size);
//  int nrepel_get_size(NoiseRepellentHandle instance);
//  void nrepel_free(NoiseRepellentHandle instance);
//  int nrepel_calculate_noise_profile(NoiseRepellentHandle
// instance,
//                                    const float *input,
//                                    float *calculated_profile);
//  int nrepel_get_noise_profile_size(NoiseRepellentHandle
// instance);  int
// nrepel_get_processing_block_size(NoiseRepellentHandle instance);
// bool nrepel_denoise(NoiseRepellentHandle instance, const float *input,
//                     float *output, const float *noise_profile,
//                     float reduction_db, float whitening_percentage,
//                     float release_ms, float noise_gain_db);
// bool nrepel_denoise_residue(NoiseRepellentHandle instance, const float
// *input,
//                             float *output, const float *noise_profile,
//                             float reduction_db, float whitening_percentage,
//                             float release_ms, float noise_gain_db);
//  bool nrepel_denoise_adaptive(NoiseRepellentHandle
// instance,
//                                    const float *input, float *output,
//                                    float reduction_db, float noise_gain_db);

#endif