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
#include <cstdint>
#include <vector>

#include <juce_audio_processors/juce_audio_processors.h>

#include "specbleach_stereo.hpp" // self-guarding for C++; do NOT wrap in extern "C"

// STFT frame-size tuning shared by the engine rebuild and the Options menu.
struct EngineConfig {
  // Stepped STFT frame-size options (ms) exposed in the Options menu.
  // Index 2 (46 ms) is the legacy default.
  static constexpr std::array<float, 5> kFrameSizeOptionsMs = {
      23.0f, 32.0f, 46.0f, 64.0f, 93.0f};
  static constexpr int kDefaultFrameSizeIndex = 2;
  // Low-latency mode: fixed 512-sample frame (~10.7 ms at 48 kHz,
  // ~11.6 ms at 44.1 kHz), causal 1D-only, zero look-ahead.
  static constexpr uint32_t kLowLatencyFrameSamples = 512;
};

// Runtime status flags shared across the audio, message and GUI threads.
// Plain struct with public atomic fields — same memory ordering and same
// thread ownership as the former NoiseRepellentAudioProcessor members.
struct ProcessorStatus {
  // Last rate reported by the host via prepareToPlay(). Zero until the first
  // call — readers must tolerate "host hasn't informed us yet" (see the
  // sr <= 0 fallback in ReductionCurveState::interpolate).
  double currentSampleRate = 0.0;
  // Consecutive silent input blocks. The engine may sleep on silence only
  // after the streak outlasts the full pipeline (latency + gain-release
  // tails) — sleeping earlier would chop decaying tails and swallow sparse
  // transients mid-response. Audio thread only.
  uint32_t silentBlocksStreak = 0;
  // Frame size (ms) the live engine group was built with. Mirrors the
  // "frame_size" APVTS choice; compared in EngineRuntime::ensureInitialized
  // to detect a rebuild request. Message thread only (written under
  // suspendProcessing).
  float currentFrameSizeMs = EngineConfig::kFrameSizeOptionsMs
                                [EngineConfig::kDefaultFrameSizeIndex];
  bool currentLowLatency = false;
  // Coalesced structural-rebuild request (frame_size/low_latency touched
  // off the message thread, e.g. LV2 automation). Consumed in processBlock.
  std::atomic<bool> frameSizeRebuildPending{false};
  // Set while setStateInformation replaces the APVTS state: a structural
  // parameter change coming from that restore must not arm the frame-size
  // profile-drop gate (GH #221 — off-thread LV2 restore defers the rebuild
  // until after the restored profile is already installed).
  std::atomic<bool> replacingStateForRestore{false};
  std::atomic<bool> preserveProfilesOnFrameRebuild{false};
  std::atomic<float> transientActivity{0.0f};
  std::atomic<bool> transientProtectionActive{false};
  std::atomic<bool> hostTransportPlaying{true};
  std::atomic<bool> hostTransportStateKnown{false};
  std::atomic<uint64_t> nonRealtimeBlockCount{0};
  uint32_t currentLatency = 0;
  uint32_t lastReportedLatency = 0; // constant per frame size / smoothing mode
  uint32_t preparedBlockSize = 0;
  uint32_t preparedNumChannels = 0;

  bool isHostTransportPlaying() const {
    return hostTransportPlaying.load();
  }
  bool isHostTransportStateKnown() const {
    return hostTransportStateKnown.load();
  }
  float consumeTransientActivity() {
    return transientActivity.exchange(0.0f);
  }
  // Counts blocks processed while the host reports non-realtime mode.
  // Some hosts (e.g. Audacity pre-initializes VST3 with kOffline at load)
  // leave isNonRealtime() stuck on while idle, so the UI gates its offline
  // indicator on this counter advancing, not on the raw flag.
  uint64_t getNonRealtimeBlockCount() const {
    return nonRealtimeBlockCount.load();
  }
};

// Thread-safe reduction-curve state: spline nodes plus the interpolated
// per-bin bias curve and the dirty flag that drives its per-block refresh.
struct ReductionCurveState {
  // Reduction Curve — spline node data
  struct CurveNode {
    float normX = 0.0f;  // 0.0 (DC) to 1.0 (Nyquist)
    float biasDB = 0.0f; // dB offset from baseline (positive = more reduction)
  };

  const std::vector<CurveNode>& nodes() const {
    return curveNodes;
  }
  void setNodes(const std::vector<CurveNode>& nodes);
  void resetNodes();
  void addNode(float normX, float biasDB);
  void removeNode(int index); // Cannot remove boundary nodes (0 or last)
  void updateNode(int index, float normX, float biasDB);

  // Callers hold lock(). The heavy lifting lives in the file-local
  // interpolateCurveInto so the setup and per-block refresh paths share it.
  void interpolate(uint32_t numBins, double sampleRate);
  // Rebuilds the interpolated reduction curve when the nodes changed. Runs
  // from the per-block parameter snapshot: the try-lock stays off the GUI's
  // node-edit path and re-arms the dirty flag when the lock is busy.
  void refreshIfDirty(bool curveEnabled, specbleach_stereo* group,
                      double sampleRate);

  const float* biasData() const {
    return interpolatedCurveBias.data();
  }
  uint32_t biasSize() const {
    return static_cast<uint32_t>(interpolatedCurveBias.size());
  }
  juce::SpinLock& lock() noexcept {
    return curveLock;
  }

  // Serializes the nodes into state's REDUCTION_CURVE child.
  void appendToState(juce::ValueTree& state);
  // Extracts REDUCTION_CURVE nodes if present (state tree first, then the raw
  // XML fallback).
  static bool decodeNodesFrom(const juce::ValueTree& state,
                              const juce::XmlElement& xmlState,
                              std::vector<CurveNode>& out);

private:
  // Thread-safe reduction curve data
  juce::SpinLock curveLock;
  std::vector<CurveNode> curveNodes{{0.0f, 0.0f}, {1.0f, 0.0f}};
  std::vector<float> interpolatedCurveBias;
  std::atomic<bool> curveNodesDirty{true};
};

// Learned noise profiles: the pending entries retained across engine
// (re)builds plus the decode/encode logic for session state.
struct ProfileStore {
  // Snapshot of one learned noise profile (per channel and mode). Shares its
  // shape with PendingProfile.
  struct SavedProfile {
    int channel = 0; // 0 = Left, 1 = Right
    int mode = 1;
    uint32_t size = 0;
    uint32_t blockCount = 0;
    std::vector<float> data;
  };

  struct PendingProfile {
    int channel = 0; // 0 = Left, 1 = Right
    int mode = 1;
    uint32_t size = 0;
    uint32_t blockCount = 0;
    std::vector<float> data;
  };

  juce::SpinLock& lock() noexcept {
    return profileLock;
  }
  std::vector<PendingProfile>& entries() noexcept {
    return pendingProfiles;
  }
  void clear() noexcept {
    pendingProfiles.clear();
  }

  bool has(specbleach_stereo* group) const;
  void reset(specbleach_stereo* group);

  // Snapshots every per-channel/per-mode profile the live group holds, so an
  // engine rebuild can restore them.
  static void captureLive(specbleach_stereo* group,
                          std::vector<SavedProfile>& out);
  // Restores a snapshot plus the retained pending entries into a rebuilt
  // group (the library resamples to the native spectrum size when the sample
  // rate changed). Entries that do not fit the current width (e.g.
  // right-channel data on a mono rebuild) are retained in pendingProfiles
  // for a later stereo rebuild.
  void restoreIntoNewGroup(specbleach_stereo* group,
                           const std::vector<SavedProfile>& savedProfiles);
  // Decodes LEARNED_PROFILES and restores it into the live group (narrow (or
  // not yet created) engines keep what they cannot take pending). Callers
  // hold lock().
  void restoreFromTree(const juce::ValueTree& profilesTree,
                       specbleach_stereo* group);
  // Serializes every channel/mode profile into state's LEARNED_PROFILES
  // child, filling per-channel/per-mode gaps before serialization.
  void appendToState(specbleach_stereo* group, juce::ValueTree& state) const;
  // Extracts LEARNED_PROFILES (state tree first, then the raw XML fallback).
  static juce::ValueTree findProfilesTree(const juce::ValueTree& state,
                                          const juce::XmlElement& xmlState);

private:
  std::vector<PendingProfile> pendingProfiles;

  // Lock to protect noise profile synchronization and serialization
  juce::SpinLock profileLock;
};

// DSP Engine — single libspecbleach multi-channel group of unified
// spectral denoisers, plus the parameter cache that gates setup-only library
// calls. The smoothing strategy (1D temporal vs 2D NLM vs 2D NLM + DFTT
// refinement) is selected per-block through parameters.smoothing_mode and
// switched seamlessly by the library (constant latency, no allocations, no
// crossfade machinery needed in the plugin).
struct EngineRuntime {
  // Parameter snapshot for one engine run (built from APVTS per block).
  struct EngineParams {
    SpecbleachDenoiserParameters p{};
    bool curveEnabled = false;
    bool learnNoise = false;
    bool adaptiveNoise = false;
    bool transientProtectionEnable = false;
  };

  EngineRuntime(ProcessorStatus& statusIn, ReductionCurveState& curveIn,
                ProfileStore& profilesIn)
      : status(statusIn), curve(curveIn), profiles(profilesIn) {
  }

  specbleach_stereo* group() const noexcept {
    return engineGroup.get();
  }
  // Transient intensity (aggregated max across channels) reported by the
  // engine for the last processed block.
  float transientIntensity() const {
    return engineGroup != nullptr
               ? specbleach_stereo_get_transient_intensity(engineGroup.get())
               : 0.0f;
  }

  EngineParams buildParams(juce::AudioProcessorValueTreeState& parameters);
  void run(juce::AudioBuffer<float>& buffer, const EngineParams& ep);
  // (Re)builds the engine group when the sample rate, channel count, frame
  // size or low-latency mode changed, restoring profiles across the rebuild.
  void ensureInitialized(double sampleRate, uint32_t wantedChannels,
                         bool lowLatency, float frameSizeMs);

private:
  // Pushes p to the engine group only when it differs from the last pushed
  // block (same-thread, serialized with process calls — never concurrent).
  // Steady-state audio blocks skip the setup-only library call entirely.
  void loadParametersIfChanged(const SpecbleachDenoiserParameters& p,
                               bool curveEnabled);

  specbleach::StereoGroupPtr engineGroup;

  // Last parameter block pushed to the engine group, with the curve pointer
  // normalized out (curve bytes are cached separately below).
  SpecbleachDenoiserParameters lastLoadedParams{};
  std::vector<float> lastLoadedCurve;
  bool paramsCacheValid = false;

  ProcessorStatus& status;
  ReductionCurveState& curve;
  ProfileStore& profiles;
};
