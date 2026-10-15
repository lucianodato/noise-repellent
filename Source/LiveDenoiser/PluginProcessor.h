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

#include <array>
#include <atomic>
#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>
#include <vector>

#include "specbleach_live_denoiser.h"

/**
 * Zero-latency time-domain noise reduction plugin.
 *
 * Wraps libspecbleach's live denoiser: a 128-band Bark-spaced time-domain
 * multiband gate. No profiles, no engine switching, no silence gating:
 * algorithmic latency is always zero, so getLatencySamples() must never
 * change after prepareToPlay.
 */
class NoiseRepellentLiveAudioProcessor
    : public juce::AudioProcessor,
      private juce::AudioProcessorValueTreeState::Listener,
      private juce::AsyncUpdater {
public:
  NoiseRepellentLiveAudioProcessor();
  ~NoiseRepellentLiveAudioProcessor() override;

  void prepareToPlay(double sampleRate, int samplesPerBlock) override;
  void releaseResources() override;

  bool isBusesLayoutSupported(const BusesLayout& layouts) const override;

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
  }
  const juce::String getProgramName(int) override {
    return {};
  }
  void changeProgramName(int, const juce::String&) override {
  }

  void getStateInformation(juce::MemoryBlock& destData) override;
  void setStateInformation(const void* data, int sizeInBytes) override;

  juce::AudioProcessorValueTreeState& getAPVTS() {
    return parameters;
  }

  // Filterbank visualization constants (matches the live engine)
  static constexpr size_t kNumBands = SPECBLEACH_LIVE_NUM_BANDS;

  // Band frame shared with GUI via lock-free ring buffer (linear
  // amplitudes straight from the engine's filterbank: input energy,
  // post-gate energy, and the gate threshold being applied).
  struct BandFrame {
    std::array<float, kNumBands> inputLevels{};
    std::array<float, kNumBands> outputLevels{};
    std::array<float, kNumBands> thresholdLevels{};
  };

  bool getNextBandFrame(BandFrame& frame);

  // Learn toggle (message thread): routes through the automatable
  // "learning" parameter. While engaged the engine tracker converges
  // on the input (floor reset on the 0->1 edge); disengaged, the
  // captured threshold freezes.
  void setLearning(bool shouldLearn);

  // Delta monitoring (message thread, GUI-only): when on, the plugin
  // outputs the removed noise (dry minus denoised) instead of the
  // denoised signal, bypassing the internal bypass.
  void setDeltaMonitoring(bool shouldMonitor);

  // Copy the engine filterbank band edges in Hz for the display scale.
  // Returns false when no engine is ready.
  bool getLiveBandEdges(float* lowerHz, float* upperHz);

  double getSampleRate() const {
    return currentSampleRate;
  }

private:
  void ensureEnginesInitialized(double sampleRate);
  void applyParameters(); // message thread only, under getCallbackLock()
  void handleAsyncUpdate() override;
  void parameterChanged(const juce::String& parameterID,
                        float newValue) override;

  juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();

  struct LiveDenoiserDeleter {
    void operator()(specbleach_live_denoiser* p) const noexcept {
      specbleach_live_denoiser_free(p);
    }
  };
  using LiveDenoiserPtr =
      std::unique_ptr<specbleach_live_denoiser, LiveDenoiserDeleter>;

  juce::AudioProcessorValueTreeState parameters;

  // One engine instance per channel (mono or stereo)
  std::array<LiveDenoiserPtr, 2> engines;
  uint32_t preparedNumChannels = 0;

  juce::AudioParameterBool* bypassParameter = nullptr;
  juce::dsp::DryWetMixer<float> dryWetMixer;
  double currentSampleRate = 44100.0;
  std::atomic<bool> parametersDirty{false};
  std::atomic<bool> learning{false};
  std::atomic<bool> deltaMonitoring{false};

  // Filterbank snapshots for GUI visualization
  juce::AbstractFifo spectralFifo{16};
  std::vector<BandFrame> spectralBuffer{16};

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentLiveAudioProcessor)
};
