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
 * Spectrum display for the live plugin, aggregated on the same ERB band
 * scale the DSAF-MP processor uses internally (38 ERB bands, 512-pt analysis
 * FFT): input (filled area), output (line), optional delta (reduction) curve,
 * and a learnable noise floor shape.
 */
class LiveSpectralVisualizerComponent : public juce::Component,
                                        public juce::Timer {
public:
  explicit LiveSpectralVisualizerComponent(
      NoiseRepellentLiveAudioProcessor& processorToUse);
  ~LiveSpectralVisualizerComponent() override;

  void paint(juce::Graphics& g) override;
  void timerCallback() override;
  void parentHierarchyChanged() override;

  // Noise floor learning: while learning, track a per-band minimum
  // statistics envelope (loop a noise-only section to converge).
  void startLearning();
  void stopLearning();
  bool isLearning() const {
    return learning;
  }

  // Delta view: show input - output (reduction applied per band)
  void setDeltaVisible(bool visible) {
    deltaVisible = visible;
  }

  static constexpr size_t kNumBands = 38; // ERB band count (matches engine)

private:
  NoiseRepellentLiveAudioProcessor& processor;
  NoiseRepellentLiveAudioProcessor::SpectralFrame currentFrame;

  std::array<float, NoiseRepellentLiveAudioProcessor::kFftBins> smoothedInputDB;
  std::array<float, NoiseRepellentLiveAudioProcessor::kFftBins>
      smoothedOutputDB;
  bool isSmoothedInitialized = false;
  int idleTicks = 0;

  bool learning = false;
  bool hasLearnedFloor = false;
  bool deltaVisible = false;

  // ERB-band aggregated values (per band, dB)
  std::array<float, kNumBands> bandInputDB{};
  std::array<float, kNumBands> bandOutputDB{};
  std::array<float, kNumBands> bandDeltaDB{};
  std::array<float, kNumBands> bandFloorDB{};
  std::array<uint32_t, kNumBands> bandStartBin{};
  std::array<uint32_t, kNumBands> bandEndBin{}; // exclusive
  uint32_t numValidBands = 0;
  uint32_t sampleRateForBands = 0;

  void rebuildBandMapping();
  void aggregateBands(
      const std::array<float, NoiseRepellentLiveAudioProcessor::kFftBins>&
          spectrumDB,
      std::array<float, kNumBands>& bandsDB) const;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LiveSpectralVisualizerComponent)
};
