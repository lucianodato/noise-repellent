/*
noise-repellent-live -- Zero-Latency Noise Reduction JUCE Plugin

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include "../PluginProcessor.h"
#include <juce_gui_basics/juce_gui_basics.h>

/**
 * RX-style band display for the live plugin. Draws the engine filterbank
 * levels directly (no FFT): input energy (filled area), post-gate output
 * energy (line), and the gate threshold being applied (the learned noise
 * floor times the aggressiveness multiplier). Learn re-arms the engine
 * noise floor so it re-converges on the current input.
 */
class LiveSpectralVisualizerComponent : public juce::Component,
                                        public juce::Timer {
public:
  explicit LiveSpectralVisualizerComponent(
      NoiseRepellentLiveAudioProcessor& processorToUse);
  ~LiveSpectralVisualizerComponent() override;

  void paint(juce::Graphics& g) override;
  void timerCallback() override;

  static constexpr size_t kNumBands =
      NoiseRepellentLiveAudioProcessor::kNumBands;

private:
  NoiseRepellentLiveAudioProcessor& processor;

  // Display-smoothed levels (dB)
  std::array<float, kNumBands> smoothedInputDB{};
  std::array<float, kNumBands> smoothedOutputDB{};
  std::array<float, kNumBands> smoothedThresholdDB{};
  bool isSmoothedInitialized = false;
  int idleTicks = 0;

  // Engine filterbank scale cache (Hz edges define the axis span)
  std::array<float, kNumBands> bandLoHz{};
  std::array<float, kNumBands> bandHiHz{};
  size_t numActiveBands = 0;
  float axisBarkLo = 0.0f;
  float axisBarkHi = 1.0f;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LiveSpectralVisualizerComponent)
};
