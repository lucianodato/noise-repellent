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

#include "VisualizationFeed.h"

#include "ProcessorSubsystems.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <span>

namespace {

using SpectralFrame = VisualizationFeed::SpectralFrame;
using VizState = VisualizationFeed::VizState;

constexpr size_t kFftSize = VisualizationFeed::kFftSize;
constexpr size_t kFftBins = VisualizationFeed::kFftBins;
constexpr size_t kMaxTonalPeaks = VisualizationFeed::kMaxTonalPeaks;

// Frame flags mirrored from the UI-facing parameters (shared by the silent and
// FFT visualization paths).
void fillFrameControlFlags(SpectralFrame& frame,
                           const juce::AudioProcessorValueTreeState& apvts) {
  frame.isLinked = (apvts.getRawParameterValue("link_reduction")->load() > 0.5f);
  frame.isOffsetLinked =
      (apvts.getRawParameterValue("link_threshold_offset")->load() > 0.5f);
  frame.reductionCurveEnabled =
      (apvts.getRawParameterValue("reduction_curve_enabled")->load() > 0.5f);
}

// Forward-transforms one accumulated window into dB magnitudes on the kFftBins
// visualization grid.
void computeSpectrumDb(const std::array<float, kFftSize>& accum,
                       std::array<float, kFftSize * 2>& work,
                       const juce::dsp::WindowingFunction<float>& window,
                       const juce::dsp::FFT& fft,
                       std::array<float, kFftBins>& magnitudeDB) {
  std::memcpy(work.data(), accum.data(), kFftSize * sizeof(float));
  std::fill(work.begin() + kFftSize, work.end(), 0.0f);
  window.multiplyWithWindowingTable(work.data(), kFftSize);
  fft.performFrequencyOnlyForwardTransform(work.data());
  for (size_t i = 0; i < kFftBins; ++i) {
    const float mag = work[i] / static_cast<float>(kFftBins);
    magnitudeDB[i] = 20.0f * std::log10(std::max(mag, 1e-7f));
  }
}

// Resamples a library noise profile (profileBins from DC 0 Hz to Nyquist
// Fs/2) onto the kFftBins visualization grid and converts to dB.
// Unnormalized power scaling offset: 20 * log10(N/2).
void fillNoiseFloorFromProfile(const float* profile, size_t profileBins,
                               std::array<float, kFftBins>& noiseFloorDB) {
  const auto maxProfileIdx = static_cast<float>(profileBins - 1);
  const auto maxFftIdx = static_cast<float>(kFftBins - 1);
  const float dbOffset =
      (maxProfileIdx > 0.0f) ? (20.0f * std::log10(maxProfileIdx)) : 0.0f;
  for (size_t i = 0; i < kFftBins; ++i) {
    const float normPos = static_cast<float>(i) / maxFftIdx;
    const float exactP = normPos * maxProfileIdx;
    const size_t p0 = std::clamp(static_cast<size_t>(exactP),
                                 static_cast<size_t>(0), profileBins - 1);
    const size_t p1 = std::min(p0 + 1, profileBins - 1);
    const float frac = exactP - static_cast<float>(p0);
    const float interpVal = (1.0f - frac) * profile[p0] + frac * profile[p1];
    const float rawDb = 10.0f * std::log10(std::max(interpVal, 1e-12f));
    noiseFloorDB[i] = rawDb - dbOffset;
  }
}

// Query libspecbleach directly for reported tonal peak frequencies computed on
// the current noise profile.
void fillTonalPeaks(SpectralFrame& frame, const specbleach_stereo* group) {
  std::array<float, kMaxTonalPeaks> peakBuf{};
  uint32_t numPeaks = 0;
  if (group != nullptr) {
    numPeaks = specbleach_stereo_get_tonal_peaks_for_channel(
        group, 0u, peakBuf.data(), static_cast<uint32_t>(peakBuf.size()));
  }
  for (uint32_t i = 0; i < numPeaks; ++i) {
    frame.tonalPeaksHz.push_back(peakBuf[i]);
  }
}

float maxTransientInWindow(const std::array<float, kFftSize>& accumTransient) {
  float maxValue = 0.0f;
  constexpr size_t centerStart = kFftSize / 4;
  constexpr size_t centerEnd = (3 * kFftSize) / 4;
  for (size_t t = centerStart; t < centerEnd; ++t) {
    if (accumTransient[t] > maxValue) {
      maxValue = accumTransient[t];
    }
  }
  return maxValue;
}

// Snapshot of the library's noise-profile visualization state (stationary or
// live adaptive). Must run under the profile try-lock.
struct VizProfile {
  const float* data = nullptr;
  uint32_t size = 0;
  bool available = false;
  bool floorEstimateAvailable = false;
};

VizProfile queryVisualizationProfile(const specbleach_stereo* group, bool isLearning,
                                     bool adaptiveNoise) {
  VizProfile vp;
  if (group == nullptr)
    return vp;
  const bool profileHasAnyMode =
      specbleach_stereo_has_any_profile_for_channel(group, 0u);
  if (isLearning) {
    vp.data = specbleach_stereo_get_noise_profile_for_channel(
        group, 0u, SPECBLEACH_PROFILE_ROLLING_MEAN);
    vp.size = specbleach_stereo_get_noise_profile_size(group);
    vp.available = (vp.data != nullptr && vp.size > 0);
    return vp;
  }
  vp.data = specbleach_stereo_get_active_noise_profile_for_channel(group, 0u);
  vp.size = specbleach_stereo_get_noise_profile_size(group);
  // The library's active profile is a double buffer allocated with the
  // group, so it can never be null while any mode profile exists.
  vp.floorEstimateAvailable =
      vp.data != nullptr && vp.size > 0 && (profileHasAnyMode || adaptiveNoise);
  vp.available = vp.floorEstimateAvailable && profileHasAnyMode;
  return vp;
}

// Fills the profile-derived part of a full spectral frame. Runs under the
// profile try-lock so library profile reads stay serialized with state
// restores (visualization try-lock in processBlock, getStateInformation).
void fillFrameProfileSection(VizState& viz, SpectralFrame& frame) {
  juce::SpinLock::ScopedTryLockType profileTryLock(viz.profileLock);
  const bool lockHeld = profileTryLock.isLocked();
  const VizProfile profile =
      lockHeld ? queryVisualizationProfile(viz.group, viz.isLearning,
                                           viz.adaptiveNoise)
               : VizProfile{};
  frame.hasNoiseProfile = profile.available;
  frame.hasNoiseFloorEstimate =
      viz.isLearning ? profile.available : profile.floorEstimateAvailable;
  fillFrameControlFlags(frame, viz.apvts);

  // Sample-accurate alignment matching the center peak of the Hann FFT
  // analysis window
  frame.transientIntensity = maxTransientInWindow(viz.accumTransient);
  frame.isTransientProtected = (frame.transientIntensity > 0.05f);
  frame.isTransientProtectionActive = viz.transientProtectionActive;
  frame.tonalPeaksHz.clear();

  if (frame.hasNoiseFloorEstimate && profile.data != nullptr) {
    fillNoiseFloorFromProfile(profile.data, profile.size, frame.noiseFloorDB);
    fillTonalPeaks(frame, viz.group);
  } else if (lockHeld) {
    // Lock was held and no profile exists — genuinely unavailable.
    // On lock contention, retain the previous frame's noise floor.
    frame.noiseFloorDB.fill(-120.0f);
  }
}

// We have a full FFT frame — compute and push to ring buffer.
void pushFullSpectralFrame(VizState& viz) {
  int start1;
  int size1;
  int start2;
  int size2;
  viz.fifo.prepareToWrite(1, start1, size1, start2, size2);

  if (size1 > 0) {
    SpectralFrame& frame = viz.frames[static_cast<size_t>(start1)];

    // ── Input FFT ──
    computeSpectrumDb(viz.accumInput, viz.inputWork, viz.window, viz.fft,
                      frame.inputMagnitudeDB);

    // ── Output FFT ──
    computeSpectrumDb(viz.accumOutput, viz.outputWork, viz.window, viz.fft,
                      frame.outputMagnitudeDB);

    // ── Noise profile from libspecbleach (stationary or live adaptive) ──
    fillFrameProfileSection(viz, frame);
  }

  if (size1 > 0)
    viz.fifo.finishedWrite(size1 + size2);
}

// Accumulates samples into the FFT windows and pushes one spectral frame per
// full window (hop = 75%). RT-safe: fixed buffers only.
void feedVisualization(VizState& viz, std::span<const float> inputSrc,
                       std::span<const float> outputSrc) {
  const size_t copySamples = inputSrc.size();
  const size_t vBufSize = viz.delayBuffer.size();
  for (size_t s = 0; s < copySamples; ++s) {
    float inSample = inputSrc[s];
    float alignedInputSample = inSample;

    if (vBufSize > 0) {
      viz.delayBuffer[viz.delayWritePos] = inSample;
      size_t delay = viz.latency % vBufSize;
      size_t readPos = (viz.delayWritePos + vBufSize - delay) % vBufSize;
      alignedInputSample = viz.delayBuffer[readPos];
      viz.delayWritePos = (viz.delayWritePos + 1) % vBufSize;
    }

    if (viz.accumCount < kFftSize) {
      viz.accumInput[viz.accumCount] = alignedInputSample;
      viz.accumOutput[viz.accumCount] = (s < outputSrc.size()) ? outputSrc[s] : 0.0f;
      viz.accumTransient[viz.accumCount] = viz.transientIntensity;
      viz.accumCount++;
    }

    if (viz.accumCount >= kFftSize) {
      pushFullSpectralFrame(viz);

      // Shift: keep last quarter for overlap (hop = 75%)
      constexpr size_t kHop = kFftSize / 4;
      std::memmove(viz.accumInput.data(), viz.accumInput.data() + kHop,
                   (kFftSize - kHop) * sizeof(float));
      std::memmove(viz.accumOutput.data(), viz.accumOutput.data() + kHop,
                   (kFftSize - kHop) * sizeof(float));
      std::memmove(viz.accumTransient.data(), viz.accumTransient.data() + kHop,
                   (kFftSize - kHop) * sizeof(float));
      viz.accumCount = kFftSize - kHop;
    }
  }
}

// Keeps the learned noise profile line frozen while input/output fall — use
// same active-profile logic as the normal FFT path so the line does not jump
// when playback starts/stops. Reads of the shared profile are guarded by the
// same try-lock as the FFT path; if the lock is busy, retain the previous
// frame's profile data.
void fillFrozenProfileSection(VizState& viz, const ProfileStore& profiles,
                              SpectralFrame& frame) {
  if (!profiles.has(viz.group)) {
    frame.noiseFloorDB.fill(-120.0f);
    frame.hasNoiseProfile = false;
    frame.hasNoiseFloorEstimate = false;
    return;
  }
  juce::SpinLock::ScopedTryLockType profileTryLock(viz.profileLock);
  if (const bool profileLockHeld = profileTryLock.isLocked();
      !profileLockHeld || viz.group == nullptr)
    return; // If the lock was busy, retain the previous frame's profile data.
  frame.hasNoiseProfile = true;
  frame.hasNoiseFloorEstimate = true;
  // Use the library's active (morphed, tonal-scaled) profile so the
  // line stays frozen and responsive even before the next audio block.
  // The active profile is a double buffer allocated together with the
  // group's profiles, so with profiles.has(group) true it is never null
  // and its size is always positive (no defensive fallback needed).
  const float* morphedPtr =
      specbleach_stereo_get_active_noise_profile_for_channel(viz.group, 0u);
  fillNoiseFloorFromProfile(
      morphedPtr, specbleach_stereo_get_noise_profile_size(viz.group),
      frame.noiseFloorDB);
  // Tonal peaks for silent when unlinked
  if ((!frame.isLinked || !frame.isOffsetLinked) && frame.hasNoiseProfile)
    fillTonalPeaks(frame, viz.group);
}
// Silent-input visualization: skip the heavy 4096 FFT but push a cheap silent
// frame so input/output fall to -120 dB while the learned noise profile stays
// frozen (only reset clears it).
void pushSilentVisualFrame(VizState& viz, const ProfileStore& profiles) {
  if ((++viz.silenceVisualCounter % 4) != 0)
    return;
  int start1;
  int size1;
  int start2;
  int size2;
  viz.fifo.prepareToWrite(1, start1, size1, start2, size2);
  if (size1 <= 0)
    return;

  SpectralFrame& frame = viz.frames[static_cast<size_t>(start1)];
  frame.inputMagnitudeDB.fill(-120.0f);
  frame.outputMagnitudeDB.fill(-120.0f);
  fillFrameControlFlags(frame, viz.apvts);
  frame.transientIntensity = 0.0f;
  frame.isTransientProtected = false;
  frame.isTransientProtectionActive = viz.transientProtectionActive;
  frame.tonalPeaksHz.clear();
  fillFrozenProfileSection(viz, profiles, frame);
  viz.fifo.finishedWrite(1);
}

} // namespace

VisualizationFeed::VizState
VisualizationFeed::makeState(const FeedInput& in) {
  return VizState{.delayBuffer = visualizerDelayBuffer,
                  .delayWritePos = visualizerDelayWritePos,
                  .accumInput = fftAccumInput,
                  .accumOutput = fftAccumOutput,
                  .accumTransient = fftAccumTransient,
                  .accumCount = fftAccumCount,
                  .silenceVisualCounter = silenceVisualCounter,
                  .fft = fftAnalyzer,
                  .window = fftWindow,
                  .inputWork = fftInputWork,
                  .outputWork = fftOutputWork,
                  .fifo = spectralFifo,
                  .frames = spectralBuffer,
                  .group = in.group,
                  .profileLock = in.profileLock,
                  .apvts = in.apvts,
                  .latency = in.latency,
                  .isLearning = in.isLearning,
                  .adaptiveNoise = in.adaptiveNoise,
                  .transientIntensity = in.transientIntensity,
                  .transientProtectionActive = in.transientProtectionActive};
}

bool VisualizationFeed::getNextFrame(SpectralFrame& frame) {
  int start1;
  int size1;
  int start2;
  int size2;
  spectralFifo.prepareToRead(1, start1, size1, start2, size2);

  if (size1 > 0) {
    frame = spectralBuffer[static_cast<size_t>(start1)];
    spectralFifo.finishedRead(1);
    return true;
  }
  return false;
}

void VisualizationFeed::feed(const FeedInput& in, std::span<const float> inputSrc,
                             std::span<const float> outputSrc) {
  VizState viz = makeState(in);
  feedVisualization(viz, inputSrc, outputSrc);
}

void VisualizationFeed::pushSilentFrame(const FeedInput& in,
                                        const ProfileStore& profiles) {
  VizState viz = makeState(in);
  pushSilentVisualFrame(viz, profiles);
}

void VisualizationFeed::resetAccumulators() {
  fftAccumInput.fill(0.0f);
  fftAccumOutput.fill(0.0f);
  fftAccumTransient.fill(0.0f);
  fftAccumCount = 0;
}

void VisualizationFeed::prepareFrames() {
  // Pre-reserve tonal peak storage in every FIFO slot so push_back on the
  // audio thread never reallocates
  for (auto& slot : spectralBuffer)
    slot.tonalPeaksHz.reserve(kMaxTonalPeaks);
}

void VisualizationFeed::resizeDelayLine(uint32_t latency) {
  visualizerDelayBuffer.assign(
      std::max<size_t>(32768, static_cast<size_t>(latency) + 8192), 0.0f);
  visualizerDelayWritePos = 0;
}
