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
#include <atomic>
#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>
#include <vector>

extern "C" {
#include "specbleach_stereo.h"
}

class NoiseRepellentAudioProcessor : public juce::AudioProcessor {
public:
  NoiseRepellentAudioProcessor();
  ~NoiseRepellentAudioProcessor() override;

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

  // Reduction Curve — spline node data
  struct CurveNode {
    float normX = 0.0f;  // 0.0 (DC) to 1.0 (Nyquist)
    float biasDB = 0.0f; // dB offset from baseline (positive = more reduction)
  };

  const std::vector<CurveNode>& getCurveNodes() const {
    return curveNodes;
  }
  void setCurveNodes(const std::vector<CurveNode>& nodes);
  void resetCurveNodes();
  void addCurveNode(float normX, float biasDB);
  void removeCurveNode(int index); // Cannot remove boundary nodes (0 or last)
  void updateCurveNode(int index, float normX, float biasDB);

  // FFT visualization constants
  static constexpr int kFftOrder = 12;               // 2^12 = 4096 point FFT
  static constexpr size_t kFftSize = 1 << kFftOrder; // 4096
  static constexpr size_t kFftBins = kFftSize / 2; // 2048 unique frequency bins

  // Spectral frame shared with GUI via lock-free ring buffer
  struct SpectralFrame {
    std::array<float, kFftBins> inputMagnitudeDB{}; // dB spectrum of input
    std::array<float, kFftBins> noiseFloorDB{}; // dB spectrum of noise profile
    std::array<float, kFftBins> outputMagnitudeDB{}; // dB spectrum of output
    std::vector<float> tonalPeaksHz{}; // Detected tonal peak frequencies in Hz
    bool hasNoiseProfile = false;
    bool isLinked = true;
    bool isOffsetLinked = true;
    bool reductionCurveEnabled = false;
    float transientIntensity = 0.0f; // Detected transient intensity [0.0, 1.0]
    bool isTransientProtected = false;
    bool isTransientProtectionActive = false;
  };

  bool getNextSpectralFrame(SpectralFrame& frame);

  void resetNoiseProfile();
  bool hasNoiseProfile() const;

  double getSampleRate() const {
    return currentSampleRate;
  }

  float consumeTransientActivity() {
    return transientActivity.exchange(0.0f, std::memory_order_relaxed);
  }

  bool isTransientProtectionActive() const {
    return transientProtectionActive.load(std::memory_order_relaxed);
  }

  // Engine-switch feedback for the UI: 1.0 == no switch in progress;
  // while switching, ramps 0..1 across warm-up + fade.
  bool isEngineSwitching() const {
    return uiSwitchProgress.load(std::memory_order_relaxed) < 1.0f;
  }

  float getEngineSwitchProgress() const {
    return uiSwitchProgress.load(std::memory_order_relaxed);
  }

private:
  void ensureEnginesInitialized(double sampleRate);
  specbleach_stereo* activeGroupFor(int algoMode) const;
  void interpolateCurve(uint32_t numBins);

  // Thread-safe reduction curve data
  juce::SpinLock curveLock;
  std::vector<CurveNode> curveNodes{{0.0f, 0.0f}, {1.0f, 0.0f}};
  std::vector<float> interpolatedCurveBias;
  std::atomic<bool> curveNodesDirty{true};

  juce::AudioProcessorValueTreeState parameters;
  static juce::AudioProcessorValueTreeState::ParameterLayout
  createParameterLayout();

  // DSP Engines — libspecbleach multi-channel groups wrapping per-channel
  // engines for both families
  specbleach_stereo* spectralGroup = nullptr; // wraps 1D per-channel engines
  specbleach_stereo* nlmGroup = nullptr;      // wraps 2D per-channel engines

  // Permanent latency alignment constants: the shorter-latency (1D) family
  // runs through this ring ONLY during engine transitions so both streams
  // share one time origin while blending. Steady-state output stays native,
  // preserving each engine's own latency and CPU profile.
  juce::AudioBuffer<float> spectralAlignmentRing;
  int alignmentWritePos = 0;
  uint32_t spectralAlignmentDelay = 0;

  // Deferred engine-switch state machine:
  // Warming: target renders in parallel until its internal buffers fill;
  //          output still comes from the source family.
  // Fading:  aligned equal-power blend toward the warmed target.
  // Slewing: after landing on the shorter-latency family, slide the
  //          alignment tap back out between delayed/direct copies.
  enum class SwitchPhase { Steady, Warming, Fading, Slewing };
  SwitchPhase switchPhase = SwitchPhase::Steady;
  int fadeFromMode = 0;               // family being blended away from
  float fadeProgress = 1.0f;          // 0..1 toward target family
  float slewProgress = 1.0f;          // 0..1 ramp-out of alignment delay
  int warmupSamplesRemaining = 0;
  static constexpr int kFadeMs = 40;
  static constexpr int kWarmupMs = 700; // >= NLM 64-frame history depth

  // GUI feedback for the ongoing switch (progress 1 == idle)
  std::atomic<float> uiSwitchProgress{1.0f};

  juce::AudioParameterBool* bypassParameter = nullptr;
  juce::dsp::DryWetMixer<float> dryWetMixer;

  double currentSampleRate = 44100.0;
  int currentAlgoMode = 1; // Track for dynamic latency updates
  std::atomic<float> transientActivity{0.0f};
  std::atomic<bool> transientProtectionActive{false};
  bool wasLearning =
      false; // Track learn mode state transition to sync profiles

  struct PendingProfile {
    int channel = 0; // 0 = Left, 1 = Right
    int mode = 1;
    uint32_t size = 0;
    uint32_t blockCount = 0;
    std::vector<float> data;
  };

  std::vector<PendingProfile> pendingProfiles;

  // Lock to protect noise profile synchronization and serialization
  juce::SpinLock profileLock;

  // Persistent dry input copy for FFT visualization (prevents RT audio thread
  // allocation)
  std::vector<float> dryInputL;
  uint32_t preparedBlockSize = 0;
  uint32_t preparedNumChannels = 0;

  // Preallocated wet scratch buffers for rendering both engine groups during
  // transitions (also serve as overflow sinks for out-of-bus channels)
  juce::AudioBuffer<float> wetScratchA; // spectralGroup wet output
  juce::AudioBuffer<float> wetScratchB; // nlmGroup wet output

  // FFT analysis for visualization
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
  std::array<float, kFftSize> fftAccumTransient{};
  size_t fftAccumCount = 0;

  // Latency-compensated delay line for input FFT visualization
  std::vector<float> visualizerDelayBuffer;
  size_t visualizerDelayWritePos = 0;
  uint32_t currentLatency = 0;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentAudioProcessor)
};
