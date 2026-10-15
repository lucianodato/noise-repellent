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

#include "ProcessorSubsystems.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <utility>

namespace {

// Divide-by-zero guards for the reduction-curve Hermite spline.
constexpr float kSplineSlopeEpsilon = 1e-5f;
constexpr float kSplineParamEpsilon = 1e-9f;
// DFTT refinement strength doubles per this many dB of reduction amount
constexpr float kDfttStrengthDoublingReductionDb = 12.0f;

using CurveNode = ReductionCurveState::CurveNode;
using SavedProfile = ProfileStore::SavedProfile;
using PendingProfile = ProfileStore::PendingProfile;

// Bounds-checked profile restore. Returns false (without touching the engine)
// when the group is null/narrow for the requested channel, so callers can
// retain the entry for a later, wider rebuild instead of dropping it.
bool loadProfileResampledChecked(specbleach_stereo* group, int channel,
                                 const float* data, uint32_t size,
                                 uint32_t blockCount,
                                 SpecbleachProfileMode mode) {
  if (group == nullptr || data == nullptr || channel < 0 || size == 0)
    return false;
  if (static_cast<uint32_t>(channel) >=
      specbleach_stereo_get_channel_count(group))
    return false;
  return specbleach_stereo_load_noise_profile_resampled_for_channel(
      group, static_cast<uint32_t>(channel), data, size, blockCount, mode);
}

// Evaluates the piecewise Hermite spline at binNormX. The node list is kept
// normX-sorted and boundary-clamped by the mutators, so the segment always
// exists; the helper walks one monotone walk instead of a loop + fallback.
size_t bracketIndex(const std::vector<CurveNode>& nodes, float x) {
  size_t idx = 0;
  while (idx + 1 < nodes.size() && nodes[idx + 1].normX <= x)
    ++idx;
  return idx;
}

bool evaluateSpline(const std::vector<CurveNode>& sortedNodes,
                    const std::vector<float>& dX, const std::vector<float>& m,
                    float binNormX, float& biasDB) {
  if (binNormX <= sortedNodes.front().normX) {
    biasDB = sortedNodes.front().biasDB;
    return true;
  }
  if (binNormX >= sortedNodes.back().normX) {
    biasDB = sortedNodes.back().biasDB;
    return true;
  }
  const size_t i = bracketIndex(sortedNodes, binNormX);
  const float dx = dX[i];
  const float t = (dx > kSplineParamEpsilon)
                      ? (binNormX - sortedNodes[i].normX) / dx
                      : 0.0f;
  const float t2 = t * t;
  const float t3 = t2 * t;

  const float h00 = 2.0f * t3 - 3.0f * t2 + 1.0f;
  const float h10 = t3 - 2.0f * t2 + t;
  const float h01 = -2.0f * t3 + 3.0f * t2;
  const float h11 = t3 - t2;

  biasDB = h00 * sortedNodes[i].biasDB + h10 * dx * m[i] +
           h01 * sortedNodes[i + 1].biasDB + h11 * dx * m[i + 1];
  return true;
}

// Slopes for the Hermite spline (endpoints and centered secant averages).
void computeSplineSlopes(size_t numNodes, const std::vector<float>& dX,
                         const std::vector<float>& dY, std::vector<float>& m) {
  m.resize(numNodes, 0.0f);
  if (numNodes <= 2) {
    const float secant =
        (dX.front() > kSplineSlopeEpsilon) ? dY.front() / dX.front() : 0.0f;
    m[0] = secant;
    m[1] = secant;
    return;
  }
  m.front() =
      (dX.front() > kSplineSlopeEpsilon) ? dY.front() / dX.front() : 0.0f;
  m.back() = (dX.back() > kSplineSlopeEpsilon) ? dY.back() / dX.back() : 0.0f;
  for (size_t i = 1; i < numNodes - 1; ++i) {
    const float secant1 =
        (dX[i - 1] > kSplineSlopeEpsilon) ? dY[i - 1] / dX[i - 1] : 0.0f;
    const float secant2 = (dX[i] > kSplineSlopeEpsilon) ? dY[i] / dX[i] : 0.0f;
    if (secant1 * secant2 <= 0.0f) {
      m[i] = 0.0f;
    } else {
      m[i] = 0.5f * (secant1 + secant2);
    }
  }
}

// Hermite-spline interpolation of the reduction-curve nodes onto numBins
// log-spaced frequency bins. Callers hold ReductionCurveState::lock(), keep
// curveNodes normX-sorted, and pass scratch buffers sized by the mutators, so
// the per-block refresh path never allocates.
void interpolateCurveInto(const std::vector<CurveNode>& curveNodes,
                          double sampleRate,
                          std::vector<float>& interpolatedCurveBias,
                          uint32_t numBins, std::vector<float>& dX,
                          std::vector<float>& dY, std::vector<float>& m) {
  if (numBins == 0)
    return;
  if (interpolatedCurveBias.size() != numBins)
    interpolatedCurveBias.resize(numBins, 0.0f);

  // setNodes/resetNodes/addNode/updateNode always keep a boundary pair, so
  // empty and single-node fills are unreachable.
  const size_t numNodes = curveNodes.size();

  double sr = sampleRate;
  if (sr <= 0.0)
    sr = 48000.0;
  const auto nyquist = static_cast<float>(sr * 0.5);

  constexpr float minFreq = 20.0f;
  constexpr float maxFreq = 20000.0f;
  const float logMin = std::log10(minFreq);
  const float logMax = std::log10(maxFreq);
  const float logRange = logMax - logMin;

  // dX/dY are pre-sized by resizeScratch for this node count (the refresh
  // path re-arms when they lag), so this path never allocates.
  jassert(dX.size() >= numNodes - 1 && dY.size() >= numNodes - 1);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    dX[i] = curveNodes[i + 1].normX - curveNodes[i].normX;
    dY[i] = curveNodes[i + 1].biasDB - curveNodes[i].biasDB;
  }

  computeSplineSlopes(numNodes, dX, dY, m);

  for (uint32_t k = 0; k < numBins; ++k) {
    float freqHz =
        (static_cast<float>(k) / static_cast<float>(numBins - 1)) * nyquist;

    float binNormX = 0.0f;
    if (freqHz <= minFreq) {
      binNormX = 0.0f;
    } else if (freqHz >= maxFreq) {
      binNormX = 1.0f;
    } else {
      binNormX = (std::log10(freqHz) - logMin) / logRange;
    }

    float biasDB = 0.0f;
    if (evaluateSpline(curveNodes, dX, m, binNormX, biasDB))
      interpolatedCurveBias[k] = biasDB;
  }
}

// Branch the smoothing strategy with plain if/else: a nested ternary chain is
// hard to extend when a new algorithm mode lands.
SpecbleachSmoothingMode smoothingModeFor(bool lowLatency, int algoMode) {
  if (lowLatency)
    return SPECBLEACH_SMOOTHING_TEMPORAL;
  if (algoMode == 3)
    return SPECBLEACH_SMOOTHING_BM3D;
  if (algoMode == 2)
    return SPECBLEACH_SMOOTHING_NLM_2D_DFTT;
  if (algoMode == 1)
    return SPECBLEACH_SMOOTHING_NLM_2D;
  return SPECBLEACH_SMOOTHING_TEMPORAL;
}

// Loads saved profiles into the group (the library resamples to the native
// spectrum size when the sample rate changed). Entries that do not fit the
// current width are retained for a later, wider rebuild instead of dropped.
// Profile snapshots share one shape with the pending-profile type, so the
// restore helpers stay templates and deduce it at the call site.
template <typename EntryT, typename OutT>
void restoreProfilesOrRetain(specbleach_stereo* group,
                             const std::vector<EntryT>& entries,
                             std::vector<OutT>& retained) {
  for (const auto& entry : entries) {
    if (loadProfileResampledChecked(group, entry.channel, entry.data.data(),
                                    entry.size, entry.blockCount,
                                    static_cast<SpecbleachProfileMode>(
                                        entry.mode)))
      continue;
    OutT pp;
    pp.channel = entry.channel;
    pp.mode = entry.mode;
    pp.size = entry.size;
    pp.blockCount = entry.blockCount;
    pp.data = entry.data;
    retained.push_back(std::move(pp));
  }
}

template <typename EntryT>
bool hasChannelOne(const std::vector<EntryT>& entries) {
  return std::any_of(entries.begin(), entries.end(),
                     [](const EntryT& entry) { return entry.channel == 1; });
}

template <typename EntryT>
void copyChannel0ProfilesToChannel1(specbleach_stereo* group,
                                    const std::vector<EntryT>& entries) {
  for (const auto& entry : entries) {
    if (entry.channel == 0) {
      loadProfileResampledChecked(group, 1, entry.data.data(), entry.size,
                                  entry.blockCount,
                                  static_cast<SpecbleachProfileMode>(
                                      entry.mode));
    }
  }
}

// Decodes one LEARNED_PROFILES child into a profile snapshot. Returns false
// for malformed entries. hasCh1 records raw channel-1 presence (legacy
// fallback) even for entries that fail to decode.
bool decodeProfileEntry(const juce::ValueTree& node, SavedProfile& out,
                        bool& hasCh1) {
  out.channel = node.getProperty("channel", 0);
  if (out.channel == 1)
    hasCh1 = true;
  out.mode = node.getProperty("mode", 0);
  out.size =
      static_cast<uint32_t>(static_cast<int>(node.getProperty("size", 0)));
  out.blockCount = static_cast<uint32_t>(
      static_cast<int>(node.getProperty("blockCount", 0)));
  const juce::String base64Data = node.getProperty("data", "");
  if (out.size == 0 || base64Data.isEmpty() || out.mode < 1 || out.mode > 4)
    return false;
  juce::MemoryBlock mb;
  if (!mb.fromBase64Encoding(base64Data) ||
      mb.getSize() < out.size * sizeof(float))
    return false;
  out.data.resize(out.size);
  std::memcpy(out.data.data(), mb.getData(), out.size * sizeof(float));
  return true;
}

void decodeProfileEntries(const juce::ValueTree& profilesTree,
                          std::vector<SavedProfile>& entries, bool& hasCh1) {
  entries.reserve(static_cast<size_t>(profilesTree.getNumChildren()));
  for (int i = 0; i < profilesTree.getNumChildren(); ++i) {
    SavedProfile entry;
    if (decodeProfileEntry(profilesTree.getChild(i), entry, hasCh1))
      entries.push_back(std::move(entry));
  }
}

} // namespace

void ReductionCurveState::setNodes(const std::vector<CurveNode>& nodes) {
  juce::SpinLock::ScopedLockType lock(curveLock);
  curveNodes = nodes;
  if (curveNodes.empty() || curveNodes.front().normX > 0.001f)
    curveNodes.insert(curveNodes.begin(), {0.0f, 0.0f});
  if (curveNodes.back().normX < 0.999f)
    curveNodes.push_back({1.0f, 0.0f});
  std::ranges::sort(curveNodes,
                    [](const CurveNode& a, const CurveNode& b) { return a.normX < b.normX; });
  resizeScratch(curveNodes);
  curveNodesDirty = true;
}

void ReductionCurveState::resetNodes() {
  juce::SpinLock::ScopedLockType lock(curveLock);
  curveNodes = {{0.0f, 0.0f}, {1.0f, 0.0f}};
  resizeScratch(curveNodes);
  curveNodesDirty = true;
}

void ReductionCurveState::addNode(float normX, float biasDB) {
  juce::SpinLock::ScopedLockType lock(curveLock);
  normX = std::clamp(normX, 0.0f, 1.0f);
  biasDB = std::clamp(biasDB, -24.0f, 24.0f);
  curveNodes.push_back({normX, biasDB});
  std::ranges::sort(curveNodes,
                    [](const CurveNode& a, const CurveNode& b) { return a.normX < b.normX; });
  resizeScratch(curveNodes);
  curveNodesDirty = true;
}

void ReductionCurveState::removeNode(int index) {
  juce::SpinLock::ScopedLockType lock(curveLock);
  if (index > 0 && index < static_cast<int>(curveNodes.size()) - 1) {
    curveNodes.erase(curveNodes.begin() + index);
    resizeScratch(curveNodes);
    curveNodesDirty = true;
  }
}

void ReductionCurveState::updateNode(int index, float normX, float biasDB) {
  juce::SpinLock::ScopedLockType lock(curveLock);
  if (index >= 0 && index < static_cast<int>(curveNodes.size())) {
    if (index == 0) {
      curveNodes[index].normX = 0.0f;
    } else if (index == static_cast<int>(curveNodes.size()) - 1) {
      curveNodes[index].normX = 1.0f;
    } else {
      curveNodes[index].normX = std::clamp(normX, 0.001f, 0.999f);
    }
    curveNodes[index].biasDB = std::clamp(biasDB, -24.0f, 24.0f);
    std::ranges::sort(curveNodes,
                      [](const CurveNode& a, const CurveNode& b) {
                        return a.normX < b.normX;
                      });
    curveNodesDirty = true;
  }
}

void ReductionCurveState::interpolate(uint32_t numBins, double sampleRate) {
  resizeScratch(curveNodes);
  interpolateCurveInto(curveNodes, sampleRate, interpolatedCurveBias, numBins,
                       dXScratch, dYScratch, curveSlopes);
}

void ReductionCurveState::refreshIfDirty(bool curveEnabled,
                                         const specbleach_stereo* group,
                                         double sampleRate) {
  if (!curveEnabled || !curveNodesDirty.exchange(false))
    return;
  juce::SpinLock::ScopedTryLockType tryLock(curveLock);
  if (!tryLock.isLocked()) {
    curveNodesDirty = true;
    return;
  }
  // Scratch is pre-sized by the mutators/setup path; the debug jassert in
  // interpolateCurveInto trips any caller that skips it.
  uint32_t profileSize = 0;
  if (group != nullptr)
    profileSize = specbleach_stereo_get_noise_profile_size(group);
  if (profileSize > 0)
    interpolateCurveInto(curveNodes, sampleRate, interpolatedCurveBias,
                         profileSize, dXScratch, dYScratch, curveSlopes);
}

void ReductionCurveState::appendToState(juce::ValueTree& state) {
  juce::ValueTree curveTree("REDUCTION_CURVE");
  {
    juce::SpinLock::ScopedLockType curveLockGuard(curveLock);
    for (const auto& node : curveNodes) {
      juce::ValueTree nodeTree("NODE");
      nodeTree.setProperty("x", static_cast<double>(node.normX), nullptr);
      nodeTree.setProperty("y", static_cast<double>(node.biasDB), nullptr);
      curveTree.appendChild(nodeTree, nullptr);
    }
  }
  // Idempotent: mirror ProfileStore::appendToState — never append twice.
  state.removeChild(state.getChildWithName("REDUCTION_CURVE"), nullptr);
  state.appendChild(curveTree, nullptr);
}

bool ReductionCurveState::decodeNodesFrom(const juce::ValueTree& state,
                                          const juce::XmlElement& xmlState,
                                          std::vector<CurveNode>& out) {
  juce::ValueTree curveTree =
      state.isValid() ? state.getChildWithName("REDUCTION_CURVE")
                      : juce::ValueTree();
  if (!curveTree.isValid()) {
    const juce::XmlElement* xmlCurve =
        xmlState.getChildByName("REDUCTION_CURVE");
    if (xmlCurve != nullptr)
      curveTree = juce::ValueTree::fromXml(*xmlCurve);
  }
  if (!curveTree.isValid() || curveTree.getNumChildren() == 0)
    return false;
  for (int i = 0; i < curveTree.getNumChildren(); ++i) {
    const juce::ValueTree nodeTree = curveTree.getChild(i);
    CurveNode cn;
    cn.normX =
        static_cast<float>(static_cast<double>(nodeTree.getProperty("x", 0.0)));
    cn.biasDB =
        static_cast<float>(static_cast<double>(nodeTree.getProperty("y", 0.0)));
    out.push_back(cn);
  }
  return true;
}

bool ProfileStore::has(const specbleach_stereo* group) const {
  if (group == nullptr)
    return false;
  return specbleach_stereo_has_any_profile_for_channel(group, 0u);
}

void ProfileStore::reset(specbleach_stereo* group) {
  specbleach_stereo_reset_profiles(group);
  pendingProfiles.clear();
}

void ProfileStore::captureLive(const specbleach_stereo* group,
                               std::vector<SavedProfile>& out) {
  if (group == nullptr)
    return;
  const uint32_t channels = specbleach_stereo_get_channel_count(group);
  const uint32_t profileSize = specbleach_stereo_get_noise_profile_size(group);
  for (uint32_t channel = 0; channel < channels; ++channel) {
    for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
         mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
      const auto profileMode = static_cast<SpecbleachProfileMode>(mode);
      if (!specbleach_stereo_profile_available_for_channel(group, channel,
                                                           profileMode))
        continue;
      const float* profile = specbleach_stereo_get_noise_profile_for_channel(
          group, channel, profileMode);
      const uint32_t blockCount =
          specbleach_stereo_get_profile_block_count_for_channel(
              group, channel, profileMode);
      if (profile != nullptr && profileSize > 0) {
        out.emplace_back(static_cast<int>(channel), mode, profileSize, blockCount,
                         std::vector<float>(profile, profile + profileSize));
      }
    }
  }
}

void ProfileStore::restoreIntoNewGroup(
    specbleach_stereo* group, const std::vector<SavedProfile>& savedProfiles) {
  std::vector<PendingProfile> retainedProfiles;
  retainedProfiles.reserve(savedProfiles.size() + pendingProfiles.size());
  restoreProfilesOrRetain(group, savedProfiles, retainedProfiles);

  // Restore state pending profiles if loaded before prepareToPlay
  const bool hasCh1 = hasChannelOne(pendingProfiles);
  restoreProfilesOrRetain(group, pendingProfiles, retainedProfiles);
  pendingProfiles = std::move(retainedProfiles);

  // Legacy states may lack channel 1 profiles — fall back from channel 0
  // (stereo groups only; mono groups have no channel 1)
  if (!hasCh1 && specbleach_stereo_get_channel_count(group) > 1) {
    copyChannel0ProfilesToChannel1(group, pendingProfiles);
  }

  // Fill any remaining per-channel/per-mode gaps
  specbleach_stereo_sync_profiles(group);
}

void ProfileStore::restoreFromTree(const juce::ValueTree& profilesTree,
                                   specbleach_stereo* group) {
  std::vector<SavedProfile> entries;
  bool hasCh1 = false;
  decodeProfileEntries(profilesTree, entries, hasCh1);

  // Narrow (or not yet created) engines cannot take this channel —
  // keep it pending for a later, wider rebuild instead of dropping
  // it.
  restoreProfilesOrRetain(group, entries, pendingProfiles);

  // If legacy profile had no Channel 1 profile, copy Channel 0 to
  // Channel 1 (stereo groups only)
  if (!hasCh1 && specbleach_stereo_get_channel_count(group) > 1) {
    copyChannel0ProfilesToChannel1(group, pendingProfiles);
  }

  // Fill any remaining per-channel/per-mode gaps
  specbleach_stereo_sync_profiles(group);
}

void ProfileStore::appendToState(specbleach_stereo* group,
                                 juce::ValueTree& state) const {
  // Fill per-channel/per-mode gaps before serialization
  specbleach_stereo_sync_profiles(group);

  juce::ValueTree profilesTree("LEARNED_PROFILES");

  const uint32_t channelCount =
      (group != nullptr) ? specbleach_stereo_get_channel_count(group) : 0;

  for (uint32_t channel = 0; channel < channelCount; ++channel) {
    for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
         mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
      const float* profile = nullptr;
      uint32_t profileSize = 0;
      uint32_t blockCount = 0;

      if (group != nullptr &&
          specbleach_stereo_profile_available_for_channel(
              group, channel, static_cast<SpecbleachProfileMode>(mode))) {
        profile = specbleach_stereo_get_noise_profile_for_channel(
            group, channel, static_cast<SpecbleachProfileMode>(mode));
        profileSize = specbleach_stereo_get_noise_profile_size(group);
        blockCount = specbleach_stereo_get_profile_block_count_for_channel(
            group, channel, static_cast<SpecbleachProfileMode>(mode));
      }

      if (profile != nullptr && profileSize > 0) {
        juce::MemoryBlock mb(profile, profileSize * sizeof(float));
        juce::ValueTree node("CHANNEL_PROFILE");
        node.setProperty("channel", static_cast<int>(channel), nullptr);
        node.setProperty("mode", mode, nullptr);
        node.setProperty("size", static_cast<int>(profileSize), nullptr);
        node.setProperty("blockCount", static_cast<int>(blockCount), nullptr);
        node.setProperty("data", mb.toBase64Encoding(), nullptr);
        profilesTree.appendChild(node, nullptr);
      }
    }
  }

  // Idempotent: a state restored earlier by replaceState can no longer carry
  // this child (the restore path scrubs it), but never append twice.
  state.removeChild(state.getChildWithName("LEARNED_PROFILES"), nullptr);
  state.appendChild(profilesTree, nullptr);
}

juce::ValueTree
ProfileStore::findProfilesTree(const juce::ValueTree& state,
                               const juce::XmlElement& xmlState) {
  juce::ValueTree profilesTree =
      state.isValid() ? state.getChildWithName("LEARNED_PROFILES")
                      : juce::ValueTree();
  if (!profilesTree.isValid()) {
    const juce::XmlElement* xmlProfiles =
        xmlState.getChildByName("LEARNED_PROFILES");
    if (xmlProfiles != nullptr)
      profilesTree = juce::ValueTree::fromXml(*xmlProfiles);
  }
  return profilesTree; } // single-line close: gcov epilogue lines are pure noise

EngineRuntime::EngineParams
EngineRuntime::buildParams(const juce::AudioProcessorValueTreeState& parameters) {
  EngineParams ep;
  // Low-latency uses a 512-sample frame: too coarse for independent
  // tonal/broadband control, so tonal follows broadband (forced link) and
  // the smoothing mode is clamped to temporal (library clamps too).
  const bool lowLatency =
      parameters.getRawParameterValue("low_latency")->load() > 0.5f;
  const auto algoMode = static_cast<int>(
      parameters.getRawParameterValue("algorithm_mode")->load());
  ep.transientProtectionEnable =
      parameters.getRawParameterValue("transient_protection_enable")->load() >
      0.5f;
  ep.learnNoise = parameters.getRawParameterValue("learn_noise")->load() > 0.5f;
  ep.adaptiveNoise =
      parameters.getRawParameterValue("adaptive_noise")->load() > 0.5f;
  const auto adaptiveMethod = static_cast<int>(
      parameters.getRawParameterValue("adaptive_method")->load());
  const float aggressiveness =
      parameters.getRawParameterValue("aggressiveness")->load();

  const bool linkReduction =
      parameters.getRawParameterValue("link_reduction")->load() > 0.5f;
  const float masterRed =
      parameters.getRawParameterValue("reduction_amount")->load();
  // Low-latency uses a 512-sample frame: too coarse for independent
  // tonal/broadband control, so tonal follows broadband (forced link).
  const float tonalRed =
      (lowLatency || linkReduction)
          ? masterRed
          : parameters.getRawParameterValue("tonal_reduction")->load();

  const float smoothing =
      parameters.getRawParameterValue("smoothing_factor")->load();

  // Masking: apply cubic curve mapping (matches original LV2 behavior)
  const float maskingRaw =
      parameters.getRawParameterValue("masking_depth")->load() / 100.0f;
  const float masking = 1.0f - std::pow(1.0f - maskingRaw, 3.0f);

  const float whitening =
      parameters.getRawParameterValue("whitening_factor")->load();

  const bool residualListen =
      parameters.getRawParameterValue("residual_listen")->load() > 0.5f;
  const float profileOffset =
      parameters.getRawParameterValue("noise_profile_offset")->load();
  const bool linkThresholdOffset =
      parameters.getRawParameterValue("link_threshold_offset")->load() > 0.5f;
  const float tonalProfileOffsetDb =
      (lowLatency || linkThresholdOffset)
          ? profileOffset
          : parameters.getRawParameterValue("tonal_noise_profile_offset")
                ->load();
  ep.curveEnabled =
      parameters.getRawParameterValue("reduction_curve_enabled")->load() > 0.5f;

  curve.refreshIfDirty(ep.curveEnabled, engineGroup.get(),
                       status.currentSampleRate);

  // Compute linear parameters for libspecbleach
  const float reductionGain = juce::Decibels::decibelsToGain(-masterRed);
  const float tonalReductionGain = juce::Decibels::decibelsToGain(-tonalRed);
  const float smoothingNorm = smoothing / 100.0f;
  const float whiteningNorm = whitening / 100.0f;
  const float profileScale = std::pow(10.0f, profileOffset / 10.0f);
  const float tonalProfileScale = std::pow(10.0f, tonalProfileOffsetDb / 10.0f);

  // Single unified parameter block. The smoothing strategy is selected via
  // smoothing_mode: the library performs the seamless runtime switch
  // internally (allocation-free crossfade, constant latency per frame size).
  ep.p.learn_noise =
      ep.learnNoise ? SPECBLEACH_LEARN_ALL : SPECBLEACH_LEARN_OFF;
  ep.p.residual_listen = residualListen;
  ep.p.reduction_gain = reductionGain;
  // Low-latency smoothing cap: the 512 frame's coarse LF bins smear gain
  // envelopes into synth-pad artifacts at long releases, so the slider's
  // top half maps to ~130 ms max instead of 500 ms.
  ep.p.smoothing_factor = lowLatency ? smoothingNorm * 0.5f : smoothingNorm;
  ep.p.smoothing_mode = smoothingModeFor(lowLatency, algoMode);
  // Single-control orchestration: in Refinement mode the DFTT quefrency
  // threshold scales with the reduction depth — deeper reduction produces
  // more musical noise, so the refinement stage is strengthed to match
  // (12 dB reduction doubles the default strength; clamped library-side).
  ep.p.dftt_strength = (!lowLatency && algoMode == 2)
                           ? (1.0f + masterRed / kDfttStrengthDoublingReductionDb)
                           : 1.0f;
  ep.p.whitening_factor = whiteningNorm;
  ep.p.adaptive_noise = ep.adaptiveNoise;
  ep.p.noise_estimation_method =
      static_cast<SpecbleachNoiseEstimationMethod>(adaptiveMethod);
  ep.p.masking_depth = masking;
  ep.p.suppression_strength = 1.0f;
  ep.p.aggressiveness = aggressiveness;
  ep.p.tonal_reduction_gain = tonalReductionGain;
  ep.p.transient_protection_enable = ep.transientProtectionEnable;
  ep.p.noise_profile_scale = profileScale;
  ep.p.reduction_curve_bias =
      ep.curveEnabled ? curve.biasData() : nullptr;
  ep.p.reduction_curve_enabled = ep.curveEnabled;
  ep.p.reduction_curve_size = ep.curveEnabled ? curve.biasSize() : 0;
  ep.p.tonal_noise_profile_scale = tonalProfileScale;
  return ep;
}

void EngineRuntime::run(juce::AudioBuffer<float>& buffer,
                        const EngineParams& ep) {
  if (engineGroup == nullptr)
    return;
  const int numChannels = buffer.getNumChannels();
  // Process 1:1 up to the buffer width (mono hosts get a 1-engine group;
  // layouts beyond stereo are rejected in isBusesLayoutSupported).
  const uint32_t groupChannels =
      std::min<uint32_t>(specbleach_stereo_get_channel_count(engineGroup.get()),
                         static_cast<uint32_t>(numChannels));
  std::array<const float*, 2> inPtrs{};
  for (uint32_t ch = 0; ch < groupChannels && ch < 2u; ++ch)
    inPtrs[ch] = buffer.getReadPointer(static_cast<int>(ch));

  std::array<float*, 2> outPtrs{};
  for (uint32_t ch = 0; ch < groupChannels && ch < 2u; ++ch)
    outPtrs[ch] = buffer.getWritePointer(static_cast<int>(ch));

  // Push the latest parameter values ahead of processing, but only when
  // they changed since the last push — steady-state blocks skip the
  // setup-only library call (same-thread handoff, never concurrent).
  loadParametersIfChanged(ep.p, ep.curveEnabled);
  specbleach_stereo_process(engineGroup.get(),
                            static_cast<uint32_t>(buffer.getNumSamples()),
                            inPtrs.data(), outPtrs.data());
}

void EngineRuntime::loadParametersIfChanged(
    const SpecbleachDenoiserParameters& p, bool curveEnabled) {
  SpecbleachDenoiserParameters norm = p;
  norm.reduction_curve_bias = nullptr; // compared via lastLoadedCurve instead
  const auto& cached = lastLoadedParams;
  // Compare values rather than object bytes because the struct may contain padding.
  const bool sameStruct =
      paramsCacheValid && norm.learn_noise == cached.learn_noise &&
      norm.residual_listen == cached.residual_listen &&
      norm.reduction_gain == cached.reduction_gain &&
      norm.smoothing_factor == cached.smoothing_factor &&
      norm.smoothing_mode == cached.smoothing_mode &&
      norm.dftt_strength == cached.dftt_strength &&
      norm.whitening_factor == cached.whitening_factor &&
      norm.adaptive_noise == cached.adaptive_noise &&
      norm.noise_estimation_method == cached.noise_estimation_method &&
      norm.masking_depth == cached.masking_depth &&
      norm.suppression_strength == cached.suppression_strength &&
      norm.aggressiveness == cached.aggressiveness &&
      norm.tonal_reduction_gain == cached.tonal_reduction_gain &&
      norm.transient_protection_enable == cached.transient_protection_enable &&
      norm.noise_profile_scale == cached.noise_profile_scale &&
      norm.reduction_curve_bias == cached.reduction_curve_bias &&
      norm.reduction_curve_enabled == cached.reduction_curve_enabled &&
      norm.reduction_curve_size == cached.reduction_curve_size &&
      norm.tonal_noise_profile_scale == cached.tonal_noise_profile_scale;
  bool sameCurve = true;
  if (curveEnabled && p.reduction_curve_bias != nullptr &&
      p.reduction_curve_size > 0) {
    sameCurve = paramsCacheValid &&
                lastLoadedCurve.size() == p.reduction_curve_size &&
                std::memcmp(lastLoadedCurve.data(), p.reduction_curve_bias,
                            p.reduction_curve_size * sizeof(float)) == 0;
  } else {
    sameCurve = paramsCacheValid && lastLoadedCurve.empty();
  }
  if (sameStruct && sameCurve)
    return; // steady state: skip the setup-only library call entirely

  if (!specbleach_stereo_load_parameters(engineGroup.get(), &p,
                                         SPECBLEACH_PARAMETERS_SIZE)) {
    return; // keep the old cache so the next block retries
  }

  lastLoadedParams = norm;
  if (curveEnabled && p.reduction_curve_bias != nullptr &&
      p.reduction_curve_size > 0)
    lastLoadedCurve.assign(p.reduction_curve_bias,
                           p.reduction_curve_bias + p.reduction_curve_size);
  else
    lastLoadedCurve.clear();
  paramsCacheValid = true;
}

void EngineRuntime::ensureInitialized(double sampleRate,
                                      uint32_t wantedChannels, bool lowLatency,
                                      float frameSizeMs) {
  // A profile-drop exemption armed by an off-thread state restore (GH #221)
  // is consumed by any engine (re)build — including the early-outs below —
  // so it can never linger to protect a later user-initiated frame switch.
  const bool preserveForStateRestore =
      status.preserveProfilesOnFrameRebuild.exchange(false);
  // Engine width follows the bus layout (mono hosts get a 1-engine group;
  // the library supports 1 channel explicitly). Rebuild on sample-rate OR
  // channel-count change — hosts can re-layout without touching the rate.
  const uint32_t currentGroupChannels =
      (engineGroup != nullptr)
          ? specbleach_stereo_get_channel_count(engineGroup.get())
          : 0;
  const float wantedFrameSizeMs =
      lowLatency && sampleRate > 0.0
          ? (static_cast<float>(EngineConfig::kLowLatencyFrameSamples) *
             1000.0f / static_cast<float>(sampleRate))
          : frameSizeMs;
  const bool wantedLowLatency = lowLatency;
  if (const bool needNewEngines =
          (engineGroup == nullptr ||
           std::abs(status.currentSampleRate - sampleRate) > 0.001 ||
           currentGroupChannels != wantedChannels ||
           std::abs(status.currentFrameSizeMs - wantedFrameSizeMs) > 0.001f ||
           status.currentLowLatency != wantedLowLatency);
      !needNewEngines)
    return;

  const bool frameSizeChanged =
      std::abs(status.currentFrameSizeMs - wantedFrameSizeMs) > 0.001f ||
      status.currentLowLatency != wantedLowLatency;

  std::vector<SavedProfile> savedProfiles;
  ProfileStore::captureLive(engineGroup.get(), savedProfiles);

  // Free existing group before allocating the replacement
  // A frame-size switch starts clean: resampling a profile captured at a
  // different resolution works poorly, so drop it and let the user re-learn
  // at native resolution. Gated on live profiles so state restores into a
  // fresh engine (empty savedProfiles) never discard session data. A rebuild
  // deferred from an off-thread state restore (GH #221) is exempt: session
  // data must survive it even when it lands in a live engine first.
  if (frameSizeChanged && !savedProfiles.empty() && !preserveForStateRestore) {
    savedProfiles.clear();
    profiles.clear();
  }
  engineGroup.reset();

  uint32_t channels = wantedChannels;

  const auto sampleRateUint = static_cast<uint32_t>(sampleRate);
  const uint32_t initFlags =
      lowLatency ? SPECBLEACH_INIT_LOW_LATENCY : 0u;
  engineGroup = specbleach::make_stereo_group(sampleRateUint, wantedFrameSizeMs,
                                              channels, initFlags);

  status.currentSampleRate = sampleRate;
  status.currentFrameSizeMs = wantedFrameSizeMs;
  status.currentLowLatency = wantedLowLatency;
  status.preparedNumChannels = channels;

  // Restore backed-up profiles into the new group (library resamples to
  // the native spectrum size when the sample rate changed). Entries that do
  // not fit the current width (e.g. right-channel data on a mono rebuild)
  // are retained in pendingProfiles for a later stereo rebuild.
  profiles.restoreIntoNewGroup(engineGroup.get(), savedProfiles);

  // This rebuild invalidates the parameter cache, and the engine is fresh.
  // Warm up the library's internal curve copy buffer here (setup context —
  // the first curve-enabled load may allocate, later same-size loads reuse
  // it) so audio blocks never trigger that allocation.
  paramsCacheValid = false;
  lastLoadedCurve.clear();
  {
    juce::SpinLock::ScopedLockType curveGuard(curve.lock());
    const uint32_t numBins =
        specbleach_stereo_get_noise_profile_size(engineGroup.get());
    if (numBins > 0) {
      curve.interpolate(numBins, status.currentSampleRate);
      lastLoadedCurve.reserve(curve.biasSize());
      if (curve.biasSize() > 0) {
        SpecbleachDenoiserParameters warm =
            specbleach_denoiser_get_default_parameters();
        warm.reduction_curve_enabled = true;
        warm.reduction_curve_bias = curve.biasData();
        warm.reduction_curve_size = curve.biasSize();
        specbleach_stereo_load_parameters(engineGroup.get(), &warm,
                                          SPECBLEACH_PARAMETERS_SIZE);
      }
    }
  }
}
