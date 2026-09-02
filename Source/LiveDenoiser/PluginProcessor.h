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

#include "specbleach_time_denoiser.h"

/**
 * Zero-latency time-domain noise reduction plugin.
 *
 * Wraps libspecbleach's DSAF-MP time-domain denoiser. No profiles, no engine
 * switching, no silence gating: algorithmic latency is always zero, so
 * getLatencySamples() must never change after prepareToPlay.
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

  // FFT visualization constants
  static constexpr int kFftOrder = 12;               // 2^12 = 4096 point FFT
  static constexpr size_t kFftSize = 1 << kFftOrder; // 4096
  static constexpr size_t kFftBins = kFftSize / 2;   // 2048 unique bins

  // Spectral frame shared with GUI via lock-free ring buffer
  struct SpectralFrame {
    std::array<float, kFftBins> inputMagnitudeDB{};  // dB spectrum of input
    std::array<float, kFftBins> outputMagnitudeDB{}; // dB spectrum of output
  };

  bool getNextSpectralFrame(SpectralFrame& frame);

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

  struct TimeDenoiserDeleter {
    void operator()(specbleach_time_denoiser* p) const noexcept {
      specbleach_time_denoiser_free(p);
    }
  };
  using TimeDenoiserPtr =
      std::unique_ptr<specbleach_time_denoiser, TimeDenoiserDeleter>;

  juce::AudioProcessorValueTreeState parameters;

  // One engine instance per channel (mono or stereo)
  std::array<TimeDenoiserPtr, 2> engines;
  uint32_t preparedNumChannels = 0;

  juce::AudioParameterBool* bypassParameter = nullptr;
  juce::AudioParameterChoice* adaptiveMethodParameter = nullptr;
  juce::dsp::DryWetMixer<float> dryWetMixer;
  double currentSampleRate = 44100.0;
  std::atomic<bool> parametersDirty{false};

  // FFT analysis for visualization (zero latency: no delay line needed)
  juce::dsp::FFT fftAnalyzer{kFftOrder};
  juce::dsp::WindowingFunction<float> fftWindow{
      kFftSize, juce::dsp::WindowingFunction<float>::hann};
  std::array<float, kFftSize * 2>
      fftInputWork{}; // real+imag interleaved for input FFT
  std::array<float, kFftSize * 2>
      fftOutputWork{}; // real+imag interleaved for output FFT

  // Lock-free SPSC Ring Buffer for GUI visualization
  juce::AbstractFifo spectralFifo{16};
  std::vector<SpectralFrame> spectralBuffer{16};

  // Accumulation buffer for FFT (collects samples across processBlock calls)
  std::array<float, kFftSize> fftAccumInput{};
  std::array<float, kFftSize> fftAccumOutput{};
  size_t fftAccumCount = 0;
  uint32_t silenceVisualCounter = 0;

  // Persistent dry input copy for FFT visualization (no RT audio-thread
  // allocation). Zero latency: dry is already aligned with wet output.
  std::vector<float> dryInputL;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentLiveAudioProcessor)
};
