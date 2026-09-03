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

#include "specbleach_stereo.hpp" // self-guarding for C++; do NOT wrap in extern "C"
#include "specbleach_version.h"

class NoiseRepellentAudioProcessor : public juce::AudioProcessor,
                                        private juce::AudioProcessorValueTreeState::Listener {
public:
  NoiseRepellentAudioProcessor();
  ~NoiseRepellentAudioProcessor() override;

  // Stepped STFT frame-size options (ms) exposed in the Options menu.
  // Index 2 (46 ms) is the legacy default.
  static constexpr float kFrameSizeOptionsMs[5] = {23.0f, 32.0f, 46.0f,
                                                  64.0f, 93.0f};
  static constexpr int kDefaultFrameSizeIndex = 2;
  float getFrameSizeMs() const;

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
  static constexpr size_t kMaxTonalPeaks = 32; // Max peaks reported per frame

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

 private:
  void ensureEnginesInitialized(double sampleRate);
  void updateLatencyReporting();
  void rebuildForFrameSizeChange();
  void parameterChanged(const juce::String& parameterID,
                        float newValue) override;
  // Parameter snapshot for one engine run (built from APVTS per block).
  struct EngineParams {
    SpecbleachDenoiserParameters p{};
    bool curveEnabled = false;
    bool learnNoise = false;
    bool adaptiveNoise = false;
    bool transientProtectionEnable = false;
  };
  EngineParams buildEngineParams();
  void runEngine(juce::AudioBuffer<float>& buffer, const EngineParams& ep);
  void interpolateCurve(uint32_t numBins);
  // Pushes p to the engine group only when it differs from the last pushed
  // block (same-thread, serialized with process calls — never concurrent).
  // Steady-state audio blocks skip the setup-only library call entirely.
  void loadParametersIfChanged(const SpecbleachDenoiserParameters& p,
                               bool curveEnabled);

  // Thread-safe reduction curve data
  juce::SpinLock curveLock;
  // Last parameter block pushed to the engine group, with the curve pointer
  // normalized out (curve bytes are cached separately below).
  SpecbleachDenoiserParameters lastLoadedParams{};
  std::vector<float> lastLoadedCurve;
  bool paramsCacheValid = false;
  std::vector<CurveNode> curveNodes{{0.0f, 0.0f}, {1.0f, 0.0f}};
  std::vector<float> interpolatedCurveBias;
  std::atomic<bool> curveNodesDirty{true};

  juce::AudioProcessorValueTreeState parameters;
  static juce::AudioProcessorValueTreeState::ParameterLayout
  createParameterLayout();

  // DSP Engine — single libspecbleach multi-channel group of unified
  // spectral denoisers. The smoothing strategy (1D temporal vs 2D NLM) is
  // selected per-block through parameters.smoothing_mode and switched
  // seamlessly by the library (constant latency, no allocations, no
  // crossfade machinery needed in the plugin).
  specbleach::StereoGroupPtr engineGroup;

  // Last rate reported by the host via prepareToPlay(). Zero until the first
  // call — readers must tolerate "host hasn't informed us yet" (see the
  // sr <= 0 fallback in interpolateCurve).
  double currentSampleRate = 0.0;
  // Consecutive silent input blocks. The engine may sleep on silence only
  // after the streak outlasts the full pipeline (latency + gain-release
  // tails) — sleeping earlier would chop decaying tails and swallow sparse
  // transients mid-response. Audio thread only.
  uint32_t silentBlocksStreak = 0;
  // Frame size (ms) the live engine group was built with. Mirrors the
  // "frame_size" APVTS choice; compared in ensureEnginesInitialized to detect
  // a rebuild request. Message thread only (written under suspendProcessing).
  float currentFrameSizeMs = kFrameSizeOptionsMs[kDefaultFrameSizeIndex];
  std::atomic<float> transientActivity{0.0f};
  std::atomic<bool> transientProtectionActive{false};

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
  uint32_t silenceVisualCounter = 0;

  // Latency-compensated delay line for input FFT visualization
  std::vector<float> visualizerDelayBuffer;
  size_t visualizerDelayWritePos = 0;
  uint32_t currentLatency = 0;
  uint32_t lastReportedLatency = 0; // constant per frame size / smoothing mode

  juce::AudioParameterBool* bypassParameter = nullptr;
  juce::dsp::DryWetMixer<float> dryWetMixer;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentAudioProcessor)
};
