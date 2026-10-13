/*
noise-repellent -- Noise Reduction JUCE Plugin

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

#include <array>
#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "ProcessorSubsystems.h"
#include "VisualizationFeed.h"

class NoiseRepellentAudioProcessor
    : public juce::AudioProcessor,
      private juce::AudioProcessorValueTreeState::Listener {
public:
  NoiseRepellentAudioProcessor();
  ~NoiseRepellentAudioProcessor() override;

  // STFT frame-size tuning, re-exported from EngineConfig (defined in
  // ProcessorSubsystems.h) for existing callers.
  static constexpr std::array<float, 5> kFrameSizeOptionsMs =
      EngineConfig::kFrameSizeOptionsMs;
  static constexpr int kDefaultFrameSizeIndex =
      EngineConfig::kDefaultFrameSizeIndex;
  static constexpr uint32_t kLowLatencyFrameSamples =
      EngineConfig::kLowLatencyFrameSamples;

  float getFrameSizeMs() const;
  bool isLowLatency() const;

  void prepareToPlay(double sampleRate, int samplesPerBlock) override;
  void releaseResources() override;

  bool isBusesLayoutSupported(const BusesLayout& layouts) const override;

  using juce::AudioProcessor::processBlock;
  using juce::AudioProcessor::processBlockBypassed;
  void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
  void processBlockBypassed(juce::AudioBuffer<float>&,
                            juce::MidiBuffer&) override;
  juce::AudioProcessorParameter* getBypassParameter() const override;

  juce::AudioProcessorEditor* createEditor() override;
  bool hasEditor() const override {
    return true;
  }

  const juce::String getName() const override {
    return JucePlugin_Name;
  }

  bool acceptsMidi() const override {
    return false;
  }
  bool producesMidi() const override {
    return false;
  }
  bool isMidiEffect() const override {
    return false;
  }
  double getTailLengthSeconds() const override {
    return 0.0;
  }

  int getNumPrograms() override {
    return 1;
  }
  int getCurrentProgram() override {
    return 0;
  }
  void setCurrentProgram(int) override {
    // Single-state plugin: exactly one program, nothing to select.
  }
  const juce::String getProgramName(int) override {
    return {};
  }
  void changeProgramName(int, const juce::String&) override {
    // Single-state plugin: no program names to change.
  }

  void getStateInformation(juce::MemoryBlock& destData) override;
  void setStateInformation(const void* data, int sizeInBytes) override;

  juce::AudioProcessorValueTreeState& getAPVTS() {
    return parameters;
  }

  // Reduction Curve — spline node data (owned by ReductionCurveState).
  using CurveNode = ReductionCurveState::CurveNode;
  // Spectral frame shared with GUI via lock-free ring buffer (owned by
  // VisualizationFeed, which also owns the FFT visualization constants).
  using SpectralFrame = VisualizationFeed::SpectralFrame;
  static constexpr int kFftOrder = VisualizationFeed::kFftOrder;
  static constexpr size_t kFftSize = VisualizationFeed::kFftSize;
  static constexpr size_t kFftBins = VisualizationFeed::kFftBins;
  static constexpr size_t kMaxTonalPeaks = VisualizationFeed::kMaxTonalPeaks;

  void resetNoiseProfile() {
    profiles.reset(engine.group());
  }
  bool hasNoiseProfile() const {
    return profiles.has(engine.group());
  }

  // Subsystems: behavior lives on these types; their methods do not count
  // against this class's method budget. See ProcessorSubsystems.h and
  // VisualizationFeed.h for the contracts and invariants.
  ProcessorStatus status;
  ReductionCurveState curve;
  ProfileStore profiles;
  EngineRuntime engine{status, curve, profiles};
  VisualizationFeed viz;

private:
  void updateLatencyReporting();
  void rebuildForFrameSizeChange(bool onAudioThread = false);
  void parameterChanged(const juce::String& parameterID,
                        float newValue) override;

  juce::AudioProcessorValueTreeState parameters;

  // Persistent dry input copy for FFT visualization (prevents RT audio thread
  // allocation)
  std::vector<float> dryInputL;

  juce::AudioParameterBool* bypassParameter = nullptr;
  // Dry-delay headroom must exceed the worst-case engine latency:
  // 93 ms at 192 kHz needs 35712 samples (latency = 2x frame).
  // setWetLatency() silently clamps past this ceiling, which
  // misaligns dry/wet and skips on bypass toggles.
  juce::dsp::DryWetMixer<float> dryWetMixer{65536};

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentAudioProcessor)
};
