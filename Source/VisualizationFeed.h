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
#include <cstddef>
#include <cstdint>
#include <vector>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "specbleach_stereo.hpp" // self-guarding for C++; do NOT wrap in extern "C"

struct ProfileStore;

// FFT analysis feed for the GUI visualization: accumulates input/output
// samples into overlapping windows, computes magnitude spectra and pushes
// SpectralFrame snapshots to the GUI through a lock-free SPSC ring buffer.
// Every method here runs on the audio thread — fixed buffers only, no
// allocations, no locks beyond the existing try-lock scopes.
struct VisualizationFeed {
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
    bool hasNoiseFloorEstimate = false;
    bool isLinked = true;
    bool isOffsetLinked = true;
    bool reductionCurveEnabled = false;
    float transientIntensity = 0.0f; // Detected transient intensity [0.0, 1.0]
    bool isTransientProtected = false;
    bool isTransientProtectionActive = false;
  };

  // Per-block inputs for the feed (everything it touches besides its own
  // buffers). References only — the feed runs on the audio thread and must
  // never allocate.
  struct FeedInput {
    specbleach_stereo* group = nullptr;
    juce::SpinLock& profileLock;
    juce::AudioProcessorValueTreeState& apvts;
    uint32_t latency = 0;
    bool isLearning = false;
    bool adaptiveNoise = false;
    float transientIntensity = 0.0f;
    bool transientProtectionActive = false;
  };

  // Everything the visualization feed touches. References only — the feed
  // runs on the audio thread and must never allocate.
  struct VizState {
    std::vector<float>& delayBuffer;
    size_t& delayWritePos;
    std::array<float, kFftSize>& accumInput;
    std::array<float, kFftSize>& accumOutput;
    std::array<float, kFftSize>& accumTransient;
    size_t& accumCount;
    uint32_t& silenceVisualCounter;
    juce::dsp::FFT& fft;
    juce::dsp::WindowingFunction<float>& window;
    std::array<float, kFftSize * 2>& inputWork;
    std::array<float, kFftSize * 2>& outputWork;
    juce::AbstractFifo& fifo;
    std::vector<SpectralFrame>& frames;
    specbleach_stereo* group = nullptr;
    juce::SpinLock& profileLock;
    juce::AudioProcessorValueTreeState& apvts;
    uint32_t latency = 0;
    bool isLearning = false;
    bool adaptiveNoise = false;
    float transientIntensity = 0.0f;
    bool transientProtectionActive = false;
  };

  bool getNextFrame(SpectralFrame& frame);
  // Accumulates samples into the FFT windows and pushes one spectral frame per
  // full window (hop = 75%). RT-safe: fixed buffers only.
  void feed(const FeedInput& in, const float* inputSrc, const float* outputSrc,
            size_t copySamples);
  // Silent-input visualization: skip the heavy 4096 FFT but push a cheap silent
  // frame so input/output fall to -120 dB while the learned noise profile stays
  // frozen (only reset clears it).
  void pushSilentFrame(const FeedInput& in, const ProfileStore& profiles);
  // Fresh STFT history: clear visualization accumulators so the first frames
  // after a pipeline (re)build don't mix pre/post-switch spectra.
  void resetAccumulators();
  void resetAccumCount() noexcept {
    fftAccumCount = 0;
  }
  // Pre-reserve tonal peak storage in every FIFO slot so push_back on the
  // audio thread never reallocates
  void prepareFrames();
  // Latency-compensated delay line for input FFT visualization
  void resizeDelayLine(uint32_t latency);

private:
  // Assembles the shared feed state from the feed's own buffers and the
  // per-block inputs.
  VizState makeState(const FeedInput& in);

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
};
