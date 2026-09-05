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

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <cmath>
#include <cstring>

namespace {

// Silence gate: amplitude floor ~-100 dB, i.e. the amplitude-domain twin of
// libspecbleach's ESTIMATOR_SILENCE_THRESHOLD (1e-10 power).
constexpr float kSilenceAmplitudeFloor = 1e-5f;
// Divide-by-zero guards for the reduction-curve Hermite spline.
constexpr float kSplineSlopeEpsilon = 1e-5f;
constexpr float kSplineParamEpsilon = 1e-9f;
// DFTT refinement strength doubles per this many dB of reduction amount
constexpr float kDfttStrengthDoublingReductionDb = 12.0f;

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

} // namespace

NoiseRepellentAudioProcessor::NoiseRepellentAudioProcessor()
    : AudioProcessor(
          BusesProperties()
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      parameters(*this, nullptr, "PARAMETERS", createParameterLayout()),
      // Dry-delay headroom must exceed the worst-case engine latency:
      // 93 ms at 192 kHz needs 35712 samples (latency = 2x frame).
      // setWetLatency() silently clamps past this ceiling, which
      // misaligns dry/wet and skips on bypass toggles.
      dryWetMixer(65536) {
  bypassParameter = dynamic_cast<juce::AudioParameterBool*>(
      parameters.getParameter("bypass"));
  parameters.addParameterListener("frame_size", this);
  parameters.addParameterListener("low_latency", this);
  juce::ignoreUnused(specbleach_get_version_string());
  DBG(specbleach_get_version_string()); // banner is "libspecbleach X.Y.Z"
}

NoiseRepellentAudioProcessor::~NoiseRepellentAudioProcessor() {
  parameters.removeParameterListener("frame_size", this);
  parameters.removeParameterListener("low_latency", this);
  releaseResources();
}

juce::AudioProcessorParameter*
NoiseRepellentAudioProcessor::getBypassParameter() const {
  return bypassParameter;
}

float NoiseRepellentAudioProcessor::getFrameSizeMs() const {
  if (auto* choice = dynamic_cast<juce::AudioParameterChoice*>(
          parameters.getParameter("frame_size"))) {
    const int index = choice->getIndex();
    if (index >= 0 && index < 5)
      return kFrameSizeOptionsMs[static_cast<size_t>(index)];
  }
  return kFrameSizeOptionsMs[kDefaultFrameSizeIndex];
}

bool NoiseRepellentAudioProcessor::isLowLatency() const {
  if (auto* p = parameters.getRawParameterValue("low_latency"))
    return p->load() > 0.5f;
  return false;
}

void NoiseRepellentAudioProcessor::parameterChanged(
    const juce::String& parameterID, float /*newValue*/) {
  if (parameterID != "frame_size" && parameterID != "low_latency")
    return;
  // APVTS listeners run on the message thread for non-automatable params
  // touched from the Options menu; state restore happens before
  // prepareToPlay so there is nothing live to rebuild yet.
  if (juce::MessageManager::getInstance()->isThisTheMessageThread()) {
    rebuildForFrameSizeChange();
  } else {
    juce::MessageManager::callAsync([this]() { rebuildForFrameSizeChange(); });
  }
}

void NoiseRepellentAudioProcessor::rebuildForFrameSizeChange() {
  if (currentSampleRate <= 0.0 || engineGroup == nullptr)
    return; // not prepared yet — next prepareToPlay picks up the new value
  // Auto-stop an in-progress Learn: a half-rolled mean must not migrate
  // across resolutions.
  if (auto* learn = parameters.getParameter("learn_noise"))
    learn->setValueNotifyingHost(0.0f);
  suspendProcessing(true);
  ensureEnginesInitialized(currentSampleRate);
  updateLatencyReporting();
  // Fresh STFT history: clear visualization accumulators so the first frames
  // after the switch don't mix pre/post-switch spectra. The rebuilt
  // pipeline is empty, so the silence streak restarts clean.
  fftAccumInput.fill(0.0f);
  fftAccumOutput.fill(0.0f);
  fftAccumTransient.fill(0.0f);
  fftAccumCount = 0;
  silentBlocksStreak = 0;
  suspendProcessing(false);
  if (auto* editor = getActiveEditor())
    editor->repaint();
}

juce::AudioProcessorValueTreeState::ParameterLayout
NoiseRepellentAudioProcessor::createParameterLayout() {
  std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

  params.push_back(std::make_unique<juce::AudioParameterChoice>(
      "algorithm_mode", "Smoothing Quality",
      juce::StringArray{"Standard (Fast & Low CPU)",
                        "Patch-Based (High Quality)",
                        "Patch-Based + Refinement (Max Quality)"},
      0));

  params.push_back(std::make_unique<juce::AudioParameterChoice>(
      "frame_size", "STFT Frame Size",
      juce::StringArray{"23 ms", "32 ms", "46 ms", "64 ms", "93 ms"},
      NoiseRepellentAudioProcessor::kDefaultFrameSizeIndex,
      juce::AudioParameterChoiceAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "low_latency", "Low Latency", false,
      juce::AudioParameterBoolAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "transient_protection_enable", "Transient Protection", false));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "learn_noise", "Learn Noise", false));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "adaptive_noise", "Adaptive Noise", false));

  params.push_back(std::make_unique<juce::AudioParameterChoice>(
      "adaptive_method", "Estimation Method",
      juce::StringArray{"SPP-MMSE (Unbiased)", "Brandt (Trimmed Mean)",
                        "Martin (Min Statistics)"},
      2));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "aggressiveness", "Aggressiveness",
      juce::NormalisableRange<float>(-1.0f, 1.0f, 0.1f), 0.5f));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "link_reduction", "Link Reduction", true,
      juce::AudioParameterBoolAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "reduction_amount", "Master Reduction",
      juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 12.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "tonal_reduction", "Tonal Reduction",
      juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 12.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "smoothing_factor", "Smoothing",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "masking_depth", "Masking Protect",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 100.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "whitening_factor", "Whitening",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "residual_listen", "Residual Listen", false));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "bypass", "Internal Bypass", false));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "show_advanced", "Show Advanced Controls", false,
      juce::AudioParameterBoolAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "show_tooltips", "Show Tooltips", true,
      juce::AudioParameterBoolAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "show_hud", "Show HUD Status", true,
      juce::AudioParameterBoolAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "noise_profile_offset", "Noise Profile Offset",
      juce::NormalisableRange<float>(-12.0f, 12.0f, 0.1f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "link_threshold_offset", "Link Threshold Offset", true,
      juce::AudioParameterBoolAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "tonal_noise_profile_offset", "Tonal Noise Profile Offset",
      juce::NormalisableRange<float>(-12.0f, 12.0f, 0.1f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "reduction_curve_enabled", "Reduction Curve", false));

  return {params.begin(), params.end()};
}

void NoiseRepellentAudioProcessor::setCurveNodes(
    const std::vector<CurveNode>& nodes) {
  juce::SpinLock::ScopedLockType lock(curveLock);
  curveNodes = nodes;
  if (curveNodes.empty() || curveNodes.front().normX > 0.001f)
    curveNodes.insert(curveNodes.begin(), {0.0f, 0.0f});
  if (curveNodes.back().normX < 0.999f)
    curveNodes.push_back({1.0f, 0.0f});
  curveNodesDirty = true;
}

void NoiseRepellentAudioProcessor::resetCurveNodes() {
  juce::SpinLock::ScopedLockType lock(curveLock);
  curveNodes = {{0.0f, 0.0f}, {1.0f, 0.0f}};
  curveNodesDirty = true;
}

void NoiseRepellentAudioProcessor::addCurveNode(float normX, float biasDB) {
  juce::SpinLock::ScopedLockType lock(curveLock);
  normX = std::clamp(normX, 0.0f, 1.0f);
  biasDB = std::clamp(biasDB, -24.0f, 24.0f);
  curveNodes.push_back({normX, biasDB});
  std::sort(
      curveNodes.begin(), curveNodes.end(),
      [](const CurveNode& a, const CurveNode& b) { return a.normX < b.normX; });
  curveNodesDirty = true;
}

void NoiseRepellentAudioProcessor::removeCurveNode(int index) {
  juce::SpinLock::ScopedLockType lock(curveLock);
  if (index > 0 && index < static_cast<int>(curveNodes.size()) - 1) {
    curveNodes.erase(curveNodes.begin() + index);
    curveNodesDirty = true;
  }
}

void NoiseRepellentAudioProcessor::updateCurveNode(int index, float normX,
                                                   float biasDB) {
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
    std::sort(curveNodes.begin(), curveNodes.end(),
              [](const CurveNode& a, const CurveNode& b) {
                return a.normX < b.normX;
              });
    curveNodesDirty = true;
  }
}

void NoiseRepellentAudioProcessor::interpolateCurve(uint32_t numBins) {
  if (numBins == 0)
    return;
  if (interpolatedCurveBias.size() != numBins)
    interpolatedCurveBias.resize(numBins, 0.0f);

  if (curveNodes.empty()) {
    std::fill(interpolatedCurveBias.begin(), interpolatedCurveBias.end(), 0.0f);
    return;
  }

  auto sortedNodes = curveNodes;
  std::sort(
      sortedNodes.begin(), sortedNodes.end(),
      [](const CurveNode& a, const CurveNode& b) { return a.normX < b.normX; });

  size_t numNodes = sortedNodes.size();
  if (numNodes == 1) {
    std::fill(interpolatedCurveBias.begin(), interpolatedCurveBias.end(),
              sortedNodes.front().biasDB);
    return;
  }

  double sr = getSampleRate();
  if (sr <= 0.0)
    sr = 48000.0;
  float nyquist = static_cast<float>(sr * 0.5);

  constexpr float minFreq = 20.0f;
  constexpr float maxFreq = 20000.0f;
  const float logMin = std::log10(minFreq);
  const float logMax = std::log10(maxFreq);
  const float logRange = logMax - logMin;

  std::vector<float> dX(numNodes - 1);
  std::vector<float> dY(numNodes - 1);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    dX[i] = sortedNodes[i + 1].normX - sortedNodes[i].normX;
    dY[i] = sortedNodes[i + 1].biasDB - sortedNodes[i].biasDB;
  }

  std::vector<float> m(numNodes, 0.0f);
  if (numNodes > 2) {
    m.front() =
        (dX.front() > kSplineSlopeEpsilon) ? dY.front() / dX.front() : 0.0f;
    m.back() = (dX.back() > kSplineSlopeEpsilon) ? dY.back() / dX.back() : 0.0f;

    for (size_t i = 1; i < numNodes - 1; ++i) {
      float secant1 =
          (dX[i - 1] > kSplineSlopeEpsilon) ? dY[i - 1] / dX[i - 1] : 0.0f;
      float secant2 = (dX[i] > kSplineSlopeEpsilon) ? dY[i] / dX[i] : 0.0f;
      if (secant1 * secant2 <= 0.0f) {
        m[i] = 0.0f;
      } else {
        m[i] = 0.5f * (secant1 + secant2);
      }
    }
  } else {
    float secant =
        (dX.front() > kSplineSlopeEpsilon) ? dY.front() / dX.front() : 0.0f;
    m[0] = secant;
    m[1] = secant;
  }

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

    if (binNormX <= sortedNodes.front().normX) {
      interpolatedCurveBias[k] = sortedNodes.front().biasDB;
      continue;
    }
    if (binNormX >= sortedNodes.back().normX) {
      interpolatedCurveBias[k] = sortedNodes.back().biasDB;
      continue;
    }

    for (size_t i = 0; i < numNodes - 1; ++i) {
      if (binNormX >= sortedNodes[i].normX &&
          binNormX <= sortedNodes[i + 1].normX) {
        float dx = dX[i];
        float t = (dx > kSplineParamEpsilon)
                      ? (binNormX - sortedNodes[i].normX) / dx
                      : 0.0f;
        float t2 = t * t;
        float t3 = t2 * t;

        float h00 = 2.0f * t3 - 3.0f * t2 + 1.0f;
        float h10 = t3 - 2.0f * t2 + t;
        float h01 = -2.0f * t3 + 3.0f * t2;
        float h11 = t3 - t2;

        float val = h00 * sortedNodes[i].biasDB + h10 * dx * m[i] +
                    h01 * sortedNodes[i + 1].biasDB + h11 * dx * m[i + 1];
        interpolatedCurveBias[k] = val;
        break;
      }
    }
  }
}

void NoiseRepellentAudioProcessor::loadParametersIfChanged(
    const SpecbleachDenoiserParameters& p, bool curveEnabled) {
  SpecbleachDenoiserParameters norm = p;
  norm.reduction_curve_bias = nullptr; // compared via lastLoadedCurve instead
  const bool sameStruct =
      paramsCacheValid &&
      std::memcmp(&norm, &lastLoadedParams, sizeof(norm)) == 0;
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
                                         SPECBLEACH_PARAMETERS_SIZE))
    return; // keep the old cache so the next block retries

  lastLoadedParams = norm;
  if (curveEnabled && p.reduction_curve_bias != nullptr &&
      p.reduction_curve_size > 0)
    lastLoadedCurve.assign(p.reduction_curve_bias,
                           p.reduction_curve_bias + p.reduction_curve_size);
  else
    lastLoadedCurve.clear();
  paramsCacheValid = true;
}

void NoiseRepellentAudioProcessor::ensureEnginesInitialized(double sampleRate) {
  // Engine width follows the bus layout (mono hosts get a 1-engine group;
  // the library supports 1 channel explicitly). Rebuild on sample-rate OR
  // channel-count change — hosts can re-layout without touching the rate.
  const uint32_t wantedChannels = static_cast<uint32_t>(
      std::max({getTotalNumInputChannels(), getTotalNumOutputChannels(), 1}));
  const uint32_t currentGroupChannels =
      (engineGroup != nullptr)
          ? specbleach_stereo_get_channel_count(engineGroup.get())
          : 0;
  const float wantedFrameSizeMs =
      isLowLatency() && sampleRate > 0.0
          ? (static_cast<float>(kLowLatencyFrameSamples) * 1000.0f /
             static_cast<float>(sampleRate))
          : getFrameSizeMs();
  const bool wantedLowLatency = isLowLatency();
  const bool needNewEngines =
      (engineGroup == nullptr ||
       std::abs(currentSampleRate - sampleRate) > 0.001 ||
       currentGroupChannels != wantedChannels ||
       std::abs(currentFrameSizeMs - wantedFrameSizeMs) > 0.001f ||
       currentLowLatency != wantedLowLatency);

  if (!needNewEngines)
    return;

  const bool frameSizeChanged =
      std::abs(currentFrameSizeMs - wantedFrameSizeMs) > 0.001f ||
      currentLowLatency != wantedLowLatency;

  struct SavedProfileData {
    int channel; // 0 = Left, 1 = Right
    int mode;
    uint32_t size;
    uint32_t blockCount;
    std::vector<float> data;
  };

  std::vector<SavedProfileData> savedProfiles;

  if (engineGroup != nullptr) {
    const uint32_t channels =
        specbleach_stereo_get_channel_count(engineGroup.get());
    const uint32_t sz =
        specbleach_stereo_get_noise_profile_size(engineGroup.get());
    for (uint32_t channel = 0; channel < channels; ++channel) {
      for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
           mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
        if (!specbleach_stereo_profile_available_for_channel(
                engineGroup.get(), channel,
                static_cast<SpecbleachProfileMode>(mode)))
          continue;
        const float* profile = specbleach_stereo_get_noise_profile_for_channel(
            engineGroup.get(), channel,
            static_cast<SpecbleachProfileMode>(mode));
        const uint32_t blockCount =
            specbleach_stereo_get_profile_block_count_for_channel(
                engineGroup.get(), channel,
                static_cast<SpecbleachProfileMode>(mode));
        if (profile != nullptr && sz > 0) {
          savedProfiles.push_back({static_cast<int>(channel), mode, sz,
                                   blockCount,
                                   std::vector<float>(profile, profile + sz)});
        }
      }
    }
  }

  // Free existing group before allocating the replacement
  // A frame-size switch starts clean: resampling a profile captured at a
  // different resolution works poorly, so drop it and let the user re-learn
  // at native resolution. Gated on live profiles so state restores into a
  // fresh engine (empty savedProfiles) never discard session data.
  if (frameSizeChanged && !savedProfiles.empty()) {
    savedProfiles.clear();
    pendingProfiles.clear();
  }
  engineGroup.reset();

  uint32_t channels = wantedChannels;

  const uint32_t sampleRateUint = static_cast<uint32_t>(sampleRate);
  const uint32_t initFlags =
      wantedLowLatency ? SPECBLEACH_INIT_LOW_LATENCY : 0u;
  engineGroup = specbleach::make_stereo_group(sampleRateUint, wantedFrameSizeMs,
                                              channels, initFlags);

  currentSampleRate = sampleRate;
  currentFrameSizeMs = wantedFrameSizeMs;
  currentLowLatency = wantedLowLatency;
  preparedNumChannels = channels;

  // Restore backed-up profiles into the new group (library resamples to
  // the native spectrum size when the sample rate changed). Entries that do
  // not fit the current width (e.g. right-channel data on a mono rebuild)
  // are retained in pendingProfiles for a later stereo rebuild.
  std::vector<PendingProfile> retainedProfiles;
  retainedProfiles.reserve(savedProfiles.size() + pendingProfiles.size());
  for (const auto& item : savedProfiles) {
    if (!loadProfileResampledChecked(
            engineGroup.get(), item.channel, item.data.data(), item.size,
            item.blockCount, static_cast<SpecbleachProfileMode>(item.mode))) {
      PendingProfile pp;
      pp.channel = item.channel;
      pp.mode = item.mode;
      pp.size = item.size;
      pp.blockCount = item.blockCount;
      pp.data = item.data;
      retainedProfiles.push_back(std::move(pp));
    }
  }

  // Restore state pending profiles if loaded before prepareToPlay
  bool hasCh1 = false;
  for (const auto& item : pendingProfiles) {
    if (item.channel == 1)
      hasCh1 = true;

    if (!loadProfileResampledChecked(
            engineGroup.get(), item.channel, item.data.data(), item.size,
            item.blockCount, static_cast<SpecbleachProfileMode>(item.mode)))
      retainedProfiles.push_back(item);
  }
  pendingProfiles = std::move(retainedProfiles);

  // Legacy states may lack channel 1 profiles — fall back from channel 0
  // (stereo groups only; mono groups have no channel 1)
  if (!hasCh1 && specbleach_stereo_get_channel_count(engineGroup.get()) > 1) {
    for (const auto& item : pendingProfiles) {
      if (item.channel == 0) {
        loadProfileResampledChecked(
            engineGroup.get(), 1, item.data.data(), item.size, item.blockCount,
            static_cast<SpecbleachProfileMode>(item.mode));
      }
    }
  }

  // Fill any remaining per-channel/per-mode gaps
  specbleach_stereo_sync_profiles(engineGroup.get());

  // This rebuild invalidates the parameter cache, and the engine is fresh.
  // Warm up the library's internal curve copy buffer here (setup context —
  // the first curve-enabled load may allocate, later same-size loads reuse
  // it) so audio blocks never trigger that allocation.
  paramsCacheValid = false;
  lastLoadedCurve.clear();
  {
    juce::SpinLock::ScopedLockType curveGuard(curveLock);
    const uint32_t numBins =
        specbleach_stereo_get_noise_profile_size(engineGroup.get());
    if (numBins > 0) {
      interpolateCurve(numBins);
      lastLoadedCurve.reserve(interpolatedCurveBias.size());
      if (!interpolatedCurveBias.empty()) {
        SpecbleachDenoiserParameters warm{};
        warm.reduction_curve_enabled = true;
        warm.reduction_curve_bias = interpolatedCurveBias.data();
        warm.reduction_curve_size =
            static_cast<uint32_t>(interpolatedCurveBias.size());
        specbleach_stereo_load_parameters(engineGroup.get(), &warm,
                                          SPECBLEACH_PARAMETERS_SIZE);
      }
    }
  }
}

void NoiseRepellentAudioProcessor::updateLatencyReporting() {
  // Latency depends on the STFT frame size (frame + NLM look-ahead) but is
  // constant across smoothing modes (temporal is padded to the NLM
  // look-ahead), so hosts only need a re-report after a frame-size switch.
  // Low-latency mode is the second exception: causal 1D-only, zero
  // look-ahead (latency = frame alone), rebuilt via the same suspend path.
  const uint32_t reportedLatency =
      engineGroup ? specbleach_stereo_get_latency(engineGroup.get()) : 0;
  lastReportedLatency = reportedLatency;
  setLatencySamples(static_cast<int>(reportedLatency));
  dryWetMixer.setWetLatency(static_cast<float>(reportedLatency));

  currentLatency = reportedLatency;
  visualizerDelayBuffer.assign(
      std::max<size_t>(32768, static_cast<size_t>(reportedLatency) + 8192),
      0.0f);
  visualizerDelayWritePos = 0;
}

void NoiseRepellentAudioProcessor::prepareToPlay(double sampleRate,
                                                 int samplesPerBlock) {
  ensureEnginesInitialized(sampleRate);

  const int bufferCapacity = std::max(samplesPerBlock, 16384);
  preparedBlockSize = static_cast<uint32_t>(bufferCapacity);
  preparedNumChannels = static_cast<uint32_t>(
      std::max({getTotalNumInputChannels(), getTotalNumOutputChannels(), 1}));

  juce::dsp::ProcessSpec spec;
  spec.sampleRate = sampleRate;
  spec.maximumBlockSize = preparedBlockSize;
  spec.numChannels = preparedNumChannels;

  dryWetMixer.prepare(spec);
  dryWetMixer.setMixingRule(juce::dsp::DryWetMixingRule::linear);
  updateLatencyReporting();

  // Pre-allocate persistent buffers to prevent audio-thread allocations
  dryInputL.resize(static_cast<size_t>(bufferCapacity), 0.0f);
  jassert(dryInputL.size() >= static_cast<size_t>(samplesPerBlock));

  // Pre-reserve tonal peak storage in every FIFO slot so push_back on the
  // audio thread never reallocates
  for (auto& slot : spectralBuffer)
    slot.tonalPeaksHz.reserve(kMaxTonalPeaks);

  // NB: pendingProfiles is owned by ensureEnginesInitialized — retained
  // (narrow-group) entries must survive prepareToPlay for a later rebuild.

  // Reset FFT accumulation
  fftAccumInput.fill(0.0f);
  fftAccumOutput.fill(0.0f);
  fftAccumTransient.fill(0.0f);
  fftAccumCount = 0;
  // Unknown pipeline state after (re)configuration: restart the silence
  // streak so a fresh silence re-proves the flush before any sleep.
  silentBlocksStreak = 0;
}

void NoiseRepellentAudioProcessor::releaseResources() {
  dryWetMixer.reset();
}

bool NoiseRepellentAudioProcessor::isBusesLayoutSupported(
    const BusesLayout& layouts) const {
  if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono() &&
      layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
    return false;

  return layouts.getMainOutputChannelSet() == layouts.getMainInputChannelSet();
}

NoiseRepellentAudioProcessor::EngineParams
NoiseRepellentAudioProcessor::buildEngineParams() {
  EngineParams ep;
  // Low-latency uses a 512-sample frame: too coarse for independent
  // tonal/broadband control, so tonal follows broadband (forced link) and
  // the smoothing mode is clamped to temporal (library clamps too).
  const bool lowLatency =
      parameters.getRawParameterValue("low_latency")->load() > 0.5f;
  const int algoMode = static_cast<int>(
      parameters.getRawParameterValue("algorithm_mode")->load());
  ep.transientProtectionEnable =
      parameters.getRawParameterValue("transient_protection_enable")->load() >
      0.5f;
  ep.learnNoise = parameters.getRawParameterValue("learn_noise")->load() > 0.5f;
  ep.adaptiveNoise =
      parameters.getRawParameterValue("adaptive_noise")->load() > 0.5f;
  const int adaptiveMethod = static_cast<int>(
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

  if (ep.curveEnabled) {
    if (curveNodesDirty.exchange(false)) {
      juce::SpinLock::ScopedTryLockType tryLock(curveLock);
      if (tryLock.isLocked()) {
        uint32_t profileSize = 0;
        if (engineGroup != nullptr)
          profileSize =
              specbleach_stereo_get_noise_profile_size(engineGroup.get());
        if (profileSize > 0)
          interpolateCurve(profileSize);
      } else {
        curveNodesDirty = true;
      }
    }
  }

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
  ep.p.smoothing_factor = smoothingNorm;
  ep.p.smoothing_mode =
      lowLatency ? SPECBLEACH_SMOOTHING_TEMPORAL
                 : (algoMode == 2)   ? SPECBLEACH_SMOOTHING_NLM_2D_DFTT
                 : (algoMode == 1) ? SPECBLEACH_SMOOTHING_NLM_2D
                                   : SPECBLEACH_SMOOTHING_TEMPORAL;
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
  ep.p.hpss_enable = ep.transientProtectionEnable;
  ep.p.noise_profile_scale = profileScale;
  ep.p.reduction_curve_bias =
      ep.curveEnabled ? interpolatedCurveBias.data() : nullptr;
  ep.p.reduction_curve_enabled = ep.curveEnabled;
  ep.p.reduction_curve_size =
      ep.curveEnabled ? static_cast<uint32_t>(interpolatedCurveBias.size())
                      : 0;
  ep.p.tonal_noise_profile_scale = tonalProfileScale;
  return ep;
}

void NoiseRepellentAudioProcessor::runEngine(juce::AudioBuffer<float>& buffer,
                                             const EngineParams& ep) {
  if (engineGroup == nullptr)
    return;
  const int numChannels = buffer.getNumChannels();
  // Process 1:1 up to the buffer width (mono hosts get a 1-engine group;
  // layouts beyond stereo are rejected in isBusesLayoutSupported).
  const uint32_t groupChannels =
      std::min<uint32_t>(specbleach_stereo_get_channel_count(engineGroup.get()),
                         static_cast<uint32_t>(numChannels));
  const float* inPtrs[2] = {nullptr, nullptr};
  for (uint32_t ch = 0; ch < groupChannels && ch < 2u; ++ch)
    inPtrs[ch] = buffer.getReadPointer(static_cast<int>(ch));

  float* outPtrs[2] = {nullptr, nullptr};
  for (uint32_t ch = 0; ch < groupChannels && ch < 2u; ++ch)
    outPtrs[ch] = buffer.getWritePointer(static_cast<int>(ch));

  // Push the latest parameter values ahead of processing, but only when
  // they changed since the last push — steady-state blocks skip the
  // setup-only library call (same-thread handoff, never concurrent).
  loadParametersIfChanged(ep.p, ep.curveEnabled);
  specbleach_stereo_process(engineGroup.get(),
                            static_cast<uint32_t>(buffer.getNumSamples()),
                            inPtrs, outPtrs);
}

void NoiseRepellentAudioProcessor::processBlock(
    juce::AudioBuffer<float>& buffer, juce::MidiBuffer&) {
  juce::ScopedNoDenormals noDenormals;
  const int numSamples = buffer.getNumSamples();
  const int numChannels = buffer.getNumChannels();

  if (numSamples == 0 || numChannels == 0)
    return;

  const bool isBypassed =
      parameters.getRawParameterValue("bypass")->load() > 0.5f;
  const EngineParams ep = buildEngineParams();
  const bool learnNoise = ep.learnNoise;
  const bool adaptiveNoise = ep.adaptiveNoise;

  // Silence detection: when input is ~-100dB and not learning, skip heavy
  // STFT+denoise (was culprit for idle > denoising and no-playback CPU).
  bool isSilent = false;
  if (!learnNoise) {
    const float thresh = kSilenceAmplitudeFloor;
    isSilent = true;
    for (int ch = 0; ch < numChannels && isSilent; ++ch) {
      const float* d = buffer.getReadPointer(ch);
      for (int i = 0; i < numSamples; ++i) {
        if (std::abs(d[i]) > thresh) {
          isSilent = false;
          break;
        }
      }
    }
    // Keep adaptive tracking alive when enabled — don't gate it as silent
    if (adaptiveNoise)
      isSilent = false;
  }

  // Save dry input copy for FFT visualization before processing
  const size_t copySamples =
      std::min(static_cast<size_t>(numSamples), dryInputL.size());
  if (numChannels >= 1 && copySamples > 0)
    std::copy_n(buffer.getReadPointer(0), copySamples, dryInputL.begin());

  // Push dry samples into JUCE DryWetMixer before in-place DSP processing
  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);

  const bool hasProfile = hasNoiseProfile();
  // Sleep the engine on silence only after the streak outlasts the whole
  // pipeline (reported latency plus gain-release tails): the pipeline
  // holds over a full latency of future output, so an output-silence check
  // cannot tell "empty" from "still primed". Bypass always runs — see below.
  if (isSilent)
    ++silentBlocksStreak;
  else
    silentBlocksStreak = 0;
  const double tailMarginSamples =
      (currentSampleRate > 0.0 ? currentSampleRate * 0.5 : 24000.0);
  const double flushSamples =
      static_cast<double>(currentLatency) + tailMarginSamples;
  const bool pipelineFlushed =
      static_cast<double>(silentBlocksStreak) *
          static_cast<double>(std::max(numSamples, 1)) >
      flushSamples;
  const bool canSkipDenoise = isSilent && pipelineFlushed && !hasProfile &&
                              !adaptiveNoise && !learnNoise;

  if (canSkipDenoise && !isBypassed) {
    // True idle: silence persisted past every tail, nothing learned or
    // tracked. The engine sleeps; cleared output keeps the silence.
    buffer.clear();
  } else {
    // The engine ALWAYS runs while bypassed. Skipping it would leave the
    // raw current input in the buffer while the mixer's wet gain is still
    // ~1 from before the toggle — the output would jump forward by a full
    // engine latency and then sweep back over the crossfade (audible skip,
    // scaling with frame size). With the engine running, both crossfade
    // endpoints carry pipeline-delayed (t - latency) content and the
    // toggle is seamless. Bonus: NLM/STFT history stays warm, so
    // un-bypassing has no re-convergence glitch.
    runEngine(buffer, ep);
  }

  // Query transient protection status and intensity (aggregated max across
  // channels)
  float reportedTransientIntensity = 0.0f;
  if (engineGroup != nullptr)
    reportedTransientIntensity =
        specbleach_stereo_get_transient_intensity(engineGroup.get());
  transientActivity.store(reportedTransientIntensity,
                          std::memory_order_relaxed);
  transientProtectionActive.store(ep.transientProtectionEnable,
                                  std::memory_order_relaxed);

  // Apply Soft Crossfade Bypass using JUCE DryWetMixer (with dry latency
  // compensation)
  dryWetMixer.setWetMixProportion(isBypassed ? 0.0f : 1.0f);
  dryWetMixer.mixWetSamples(audioBlock);

  // Skip FFT analysis and FIFO writes during offline rendering or if the GUI
  // is closed
  if (isNonRealtime() || getActiveEditor() == nullptr) {
    fftAccumCount = 0;
    return;
  }

  // When input is silent and not learning, skip heavy 4096 FFT but push a
  // cheap silent frame so input/output fall to -120 dB while the learned
  // noise profile stays frozen (only reset clears it).
  if (isSilent) {
    fftAccumCount = 0;
    if ((++silenceVisualCounter % 4) == 0) {
      int start1, size1, start2, size2;
      spectralFifo.prepareToWrite(1, start1, size1, start2, size2);
      if (size1 > 0) {
        SpectralFrame& frame = spectralBuffer[static_cast<size_t>(start1)];
        frame.inputMagnitudeDB.fill(-120.0f);
        frame.outputMagnitudeDB.fill(-120.0f);
        frame.isLinked =
            (parameters.getRawParameterValue("link_reduction")->load() > 0.5f);
        frame.isOffsetLinked =
            (parameters.getRawParameterValue("link_threshold_offset")->load() >
             0.5f);
        frame.reductionCurveEnabled =
            (parameters.getRawParameterValue("reduction_curve_enabled")
                 ->load() > 0.5f);
        frame.transientIntensity = 0.0f;
        frame.isTransientProtected = false;
        frame.isTransientProtectionActive = ep.transientProtectionEnable;
        frame.tonalPeaksHz.clear();
        // Keep learned noise profile frozen while input/output fall — use
        // same active-profile logic as the normal FFT path so the line does
        // not jump when playback starts/stops. Reads of the shared profile
        // are guarded by the same try-lock as the FFT path; if the lock is
        // busy, retain the previous frame's profile data.
        if (hasNoiseProfile()) {
          juce::SpinLock::ScopedTryLockType profileTryLock(profileLock);
          if (profileTryLock.isLocked() && engineGroup != nullptr) {
            frame.hasNoiseProfile = true;
            // Use the library's active (morphed, tonal-scaled) profile so the
            // line stays frozen and responsive even before the next audio
            // block.
            uint32_t profileSize = 0;
            const float* morphedPtr =
                specbleach_stereo_get_active_noise_profile_for_channel(
                    engineGroup.get(), 0u);
            profileSize =
                specbleach_stereo_get_noise_profile_size(engineGroup.get());
            if (morphedPtr != nullptr && profileSize > 0) {
              // Resample morphed profile to kFftBins and convert to dB
              size_t realProfileBins = profileSize;
              const float maxProfileIdx =
                  static_cast<float>(realProfileBins - 1);
              const float maxFftIdx = static_cast<float>(kFftBins - 1);
              const float dbOffset = (maxProfileIdx > 0.0f)
                                         ? (20.0f * std::log10(maxProfileIdx))
                                         : 0.0f;
              for (size_t i = 0; i < kFftBins; ++i) {
                float normPos = static_cast<float>(i) / maxFftIdx;
                float exactP = normPos * maxProfileIdx;
                size_t p0 =
                    std::clamp(static_cast<size_t>(exactP),
                               static_cast<size_t>(0), realProfileBins - 1);
                size_t p1 = std::min(p0 + 1, realProfileBins - 1);
                float frac = exactP - static_cast<float>(p0);
                float interpVal =
                    (1.0f - frac) * morphedPtr[p0] + frac * morphedPtr[p1];
                float rawDb = 10.0f * std::log10(std::max(interpVal, 1e-12f));
                frame.noiseFloorDB[i] = rawDb - dbOffset;
              }
            } else {
              frame.noiseFloorDB.fill(-120.0f);
            }
            // Tonal peaks for silent when unlinked
            if ((!frame.isLinked || !frame.isOffsetLinked) &&
                frame.hasNoiseProfile) {
              std::array<float, kMaxTonalPeaks> peakBuf{};
              uint32_t n = specbleach_stereo_get_tonal_peaks_for_channel(
                  engineGroup.get(), 0u, peakBuf.data(),
                  static_cast<uint32_t>(peakBuf.size()));
              for (uint32_t i = 0; i < n; ++i)
                frame.tonalPeaksHz.push_back(peakBuf[i]);
            }
          }
          // If the lock was busy, retain the previous frame's profile data.
        } else {
          frame.noiseFloorDB.fill(-120.0f);
          frame.hasNoiseProfile = false;
        }
        spectralFifo.finishedWrite(1);
      }
    }
    return;
  }

  // ── FFT-based spectral frame for GUI visualization ──
  // Accumulate samples until we have a full FFT window (kFftSize samples)
  const float* inputSrc = dryInputL.data();
  const float* outputSrc =
      (numChannels >= 1) ? buffer.getReadPointer(0) : nullptr;

  const size_t vBufSize = visualizerDelayBuffer.size();
  for (size_t s = 0; s < copySamples; ++s) {
    float inSample = inputSrc[s];
    float alignedInputSample = inSample;

    if (vBufSize > 0) {
      visualizerDelayBuffer[visualizerDelayWritePos] = inSample;
      size_t delay = currentLatency % vBufSize;
      size_t readPos = (visualizerDelayWritePos + vBufSize - delay) % vBufSize;
      alignedInputSample = visualizerDelayBuffer[readPos];
      visualizerDelayWritePos = (visualizerDelayWritePos + 1) % vBufSize;
    }

    const float currentTransientVal = reportedTransientIntensity;
    if (fftAccumCount < kFftSize) {
      fftAccumInput[fftAccumCount] = alignedInputSample;
      fftAccumOutput[fftAccumCount] = outputSrc ? outputSrc[s] : 0.0f;
      fftAccumTransient[fftAccumCount] = currentTransientVal;
      fftAccumCount++;
    }

    if (fftAccumCount >= kFftSize) {
      // We have a full FFT frame — compute and push to ring buffer
      int start1, size1, start2, size2;
      spectralFifo.prepareToWrite(1, start1, size1, start2, size2);

      if (size1 > 0) {
        SpectralFrame& frame = spectralBuffer[static_cast<size_t>(start1)];

        // ── Input FFT ──
        std::memcpy(fftInputWork.data(), fftAccumInput.data(),
                    kFftSize * sizeof(float));
        std::fill(fftInputWork.begin() + kFftSize, fftInputWork.end(), 0.0f);
        fftWindow.multiplyWithWindowingTable(fftInputWork.data(), kFftSize);
        fftAnalyzer.performFrequencyOnlyForwardTransform(fftInputWork.data());

        for (size_t i = 0; i < kFftBins; ++i) {
          float mag = fftInputWork[i] / static_cast<float>(kFftBins);
          frame.inputMagnitudeDB[i] = 20.0f * std::log10(std::max(mag, 1e-7f));
        }

        // ── Output FFT ──
        std::memcpy(fftOutputWork.data(), fftAccumOutput.data(),
                    kFftSize * sizeof(float));
        std::fill(fftOutputWork.begin() + kFftSize, fftOutputWork.end(), 0.0f);
        fftWindow.multiplyWithWindowingTable(fftOutputWork.data(), kFftSize);
        fftAnalyzer.performFrequencyOnlyForwardTransform(fftOutputWork.data());

        for (size_t i = 0; i < kFftBins; ++i) {
          float mag = fftOutputWork[i] / static_cast<float>(kFftBins);
          frame.outputMagnitudeDB[i] = 20.0f * std::log10(std::max(mag, 1e-7f));
        }

        // ── Noise profile from libspecbleach (stationary or live adaptive) ──
        const float* actualNoiseProfile = nullptr;
        uint32_t profileSize = 0;
        bool profileAvailable = false;

        const bool isLearning = learnNoise;
        bool profileHasAnyMode = false;

        juce::SpinLock::ScopedTryLockType profileTryLock(profileLock);
        if (profileTryLock.isLocked() && engineGroup != nullptr) {
          const auto* vizGroup = engineGroup.get();
          profileHasAnyMode =
              specbleach_stereo_has_any_profile_for_channel(vizGroup, 0u);

          if (isLearning) {
            actualNoiseProfile =
                specbleach_stereo_get_noise_profile_for_channel(
                    vizGroup, 0u, SPECBLEACH_PROFILE_ROLLING_MEAN);
            profileSize = specbleach_stereo_get_noise_profile_size(vizGroup);
            profileAvailable =
                (actualNoiseProfile != nullptr && profileSize > 0);
          } else {
            actualNoiseProfile =
                specbleach_stereo_get_active_noise_profile_for_channel(vizGroup,
                                                                       0u);
            profileSize = specbleach_stereo_get_noise_profile_size(vizGroup);
            if (!actualNoiseProfile && profileHasAnyMode) {
              for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
                   mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
                if (specbleach_stereo_profile_available_for_channel(
                        vizGroup, 0u,
                        static_cast<SpecbleachProfileMode>(mode))) {
                  actualNoiseProfile =
                      specbleach_stereo_get_noise_profile_for_channel(
                          vizGroup, 0u,
                          static_cast<SpecbleachProfileMode>(mode));
                  if (actualNoiseProfile)
                    break;
                }
              }
            }
            profileAvailable =
                (actualNoiseProfile != nullptr && profileSize > 0);
          }
        }

        frame.hasNoiseProfile = profileAvailable;
        frame.isLinked =
            (parameters.getRawParameterValue("link_reduction")->load() > 0.5f);
        frame.isOffsetLinked =
            (parameters.getRawParameterValue("link_threshold_offset")->load() >
             0.5f);
        frame.reductionCurveEnabled =
            (parameters.getRawParameterValue("reduction_curve_enabled")
                 ->load() > 0.5f);

        // Sample-accurate alignment matching the center peak of the Hann FFT
        // analysis window
        float maxWindowTransient = 0.0f;
        constexpr size_t centerStart = kFftSize / 4;
        constexpr size_t centerEnd = (3 * kFftSize) / 4;
        for (size_t t = centerStart; t < centerEnd; ++t) {
          if (fftAccumTransient[t] > maxWindowTransient) {
            maxWindowTransient = fftAccumTransient[t];
          }
        }
        frame.transientIntensity = maxWindowTransient;
        frame.isTransientProtected = (maxWindowTransient > 0.05f);
        frame.isTransientProtectionActive = ep.transientProtectionEnable;
        frame.tonalPeaksHz.clear();

        if (profileAvailable && actualNoiseProfile) {
          // profileSize represents the unique spectrum from DC (0 Hz) to
          // Nyquist (Fs/2)
          size_t realProfileBins = profileSize;

          const float maxProfileIdx = static_cast<float>(realProfileBins - 1);
          const float maxFftIdx = static_cast<float>(kFftBins - 1);
          // Unnormalized power scaling offset: 20 * log10(N/2)
          const float dbOffset = (maxProfileIdx > 0.0f)
                                     ? (20.0f * std::log10(maxProfileIdx))
                                     : 0.0f;

          for (size_t i = 0; i < kFftBins; ++i) {
            float normPos = static_cast<float>(i) / maxFftIdx;
            float exactP = normPos * maxProfileIdx;

            size_t p0 = std::clamp(static_cast<size_t>(exactP),
                                   static_cast<size_t>(0), realProfileBins - 1);
            size_t p1 = std::min(p0 + 1, realProfileBins - 1);
            float frac = exactP - static_cast<float>(p0);

            float val0 = actualNoiseProfile[p0];
            float val1 = actualNoiseProfile[p1];
            float interpVal = (1.0f - frac) * val0 + frac * val1;

            float rawDb = 10.0f * std::log10(std::max(interpVal, 1e-12f));
            frame.noiseFloorDB[i] = rawDb - dbOffset;
          }

          // Query libspecbleach directly for reported tonal peak frequencies
          // computed on this real noise profile
          std::array<float, kMaxTonalPeaks> peakBuf{};
          uint32_t numPeaks = 0;
          if (engineGroup != nullptr) {
            numPeaks = specbleach_stereo_get_tonal_peaks_for_channel(
                engineGroup.get(), 0u, peakBuf.data(),
                static_cast<uint32_t>(peakBuf.size()));
          }

          for (uint32_t i = 0; i < numPeaks; ++i) {
            frame.tonalPeaksHz.push_back(peakBuf[i]);
          }
        } else if (profileTryLock.isLocked()) {
          // Lock was held and no profile exists — genuinely unavailable.
          // On lock contention, retain the previous frame's noise floor.
          frame.noiseFloorDB.fill(-120.0f);
        }
      }

      spectralFifo.finishedWrite(1);

      // Shift: keep last quarter for overlap (hop = 75%)
      constexpr size_t kHop = kFftSize / 4;
      std::memmove(fftAccumInput.data(), fftAccumInput.data() + kHop,
                   (kFftSize - kHop) * sizeof(float));
      std::memmove(fftAccumOutput.data(), fftAccumOutput.data() + kHop,
                   (kFftSize - kHop) * sizeof(float));
      std::memmove(fftAccumTransient.data(), fftAccumTransient.data() + kHop,
                   (kFftSize - kHop) * sizeof(float));
      fftAccumCount = kFftSize - kHop;
    }
  }
}

void NoiseRepellentAudioProcessor::processBlockBypassed(
    juce::AudioBuffer<float>& buffer, juce::MidiBuffer& /*midiMessages*/) {
  const int numSamples = buffer.getNumSamples();
  const int numChannels = buffer.getNumChannels();

  if (numSamples == 0 || numChannels == 0)
    return;

  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);
  // Run the engine even though the wet path is muted: same alignment
  // requirement as the internal bypass (see processBlock) — the buffer
  // must hold pipeline-delayed content while the mixer ramps down, and
  // the delay line stays warm for the resume boundary.
  runEngine(buffer, buildEngineParams());
  dryWetMixer.setWetMixProportion(0.0f);
  dryWetMixer.mixWetSamples(audioBlock);
}

bool NoiseRepellentAudioProcessor::getNextSpectralFrame(SpectralFrame& frame) {
  int start1, size1, start2, size2;
  spectralFifo.prepareToRead(1, start1, size1, start2, size2);

  if (size1 > 0) {
    frame = spectralBuffer[static_cast<size_t>(start1)];
    spectralFifo.finishedRead(1);
    return true;
  }
  return false;
}

void NoiseRepellentAudioProcessor::resetNoiseProfile() {
  specbleach_stereo_reset_profiles(engineGroup.get());
  pendingProfiles.clear();
}

bool NoiseRepellentAudioProcessor::hasNoiseProfile() const {
  if (engineGroup == nullptr)
    return false;
  return specbleach_stereo_has_any_profile_for_channel(engineGroup.get(), 0u);
}

void NoiseRepellentAudioProcessor::getStateInformation(
    juce::MemoryBlock& destData) {
  juce::SpinLock::ScopedLockType lock(profileLock);

  // Fill per-channel/per-mode gaps before serialization
  specbleach_stereo_sync_profiles(engineGroup.get());

  auto state = parameters.copyState();

  juce::ValueTree profilesTree("LEARNED_PROFILES");

  const uint32_t channelCount =
      (engineGroup != nullptr)
          ? specbleach_stereo_get_channel_count(engineGroup.get())
          : 0;

  for (uint32_t channel = 0; channel < channelCount; ++channel) {
    for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
         mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
      const float* profile = nullptr;
      uint32_t profileSize = 0;
      uint32_t blockCount = 0;

      if (engineGroup != nullptr &&
          specbleach_stereo_profile_available_for_channel(
              engineGroup.get(), channel,
              static_cast<SpecbleachProfileMode>(mode))) {
        profile = specbleach_stereo_get_noise_profile_for_channel(
            engineGroup.get(), channel,
            static_cast<SpecbleachProfileMode>(mode));
        profileSize =
            specbleach_stereo_get_noise_profile_size(engineGroup.get());
        blockCount = specbleach_stereo_get_profile_block_count_for_channel(
            engineGroup.get(), channel,
            static_cast<SpecbleachProfileMode>(mode));
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

  state.appendChild(profilesTree, nullptr);

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
  state.appendChild(curveTree, nullptr);

  std::unique_ptr<juce::XmlElement> xml(state.createXml());
  copyXmlToBinary(*xml, destData);
}

void NoiseRepellentAudioProcessor::setStateInformation(const void* data,
                                                       int sizeInBytes) {
  pendingProfiles.clear();

  std::unique_ptr<juce::XmlElement> xmlState(
      getXmlFromBinary(data, sizeInBytes));
  if (xmlState == nullptr && data != nullptr && sizeInBytes > 0) {
    juce::String xmlString =
        juce::String::createStringFromData(data, sizeInBytes);
    if (xmlString.trim().startsWith("<")) {
      xmlState = juce::XmlDocument::parse(xmlString);
    }
  }

  if (xmlState != nullptr) {
    juce::ValueTree state = juce::ValueTree::fromXml(*xmlState);

    // Extract REDUCTION_CURVE if present
    juce::ValueTree curveTree;
    if (state.isValid()) {
      curveTree = state.getChildWithName("REDUCTION_CURVE");
    }
    if (!curveTree.isValid()) {
      juce::XmlElement* xmlCurve = xmlState->getChildByName("REDUCTION_CURVE");
      if (xmlCurve != nullptr) {
        curveTree = juce::ValueTree::fromXml(*xmlCurve);
      }
    }
    if (curveTree.isValid() && curveTree.getNumChildren() > 0) {
      std::vector<CurveNode> loadedNodes;
      for (int i = 0; i < curveTree.getNumChildren(); ++i) {
        juce::ValueTree nodeTree = curveTree.getChild(i);
        CurveNode cn;
        cn.normX = static_cast<float>(
            static_cast<double>(nodeTree.getProperty("x", 0.0)));
        cn.biasDB = static_cast<float>(
            static_cast<double>(nodeTree.getProperty("y", 0.0)));
        loadedNodes.push_back(cn);
      }
      setCurveNodes(loadedNodes);
    }

    // Extract LEARNED_PROFILES before replaceState modifies state
    juce::ValueTree profilesTree;
    if (state.isValid()) {
      profilesTree = state.getChildWithName("LEARNED_PROFILES");
    }

    if (!profilesTree.isValid()) {
      juce::XmlElement* xmlProfiles =
          xmlState->getChildByName("LEARNED_PROFILES");
      if (xmlProfiles != nullptr) {
        profilesTree = juce::ValueTree::fromXml(*xmlProfiles);
      }
    }

    if (state.isValid()) {
      parameters.replaceState(state);
    }

    if (profilesTree.isValid()) {
      // Serialize the whole restore sequence against profile readers
      // (visualization try-lock in processBlock, getStateInformation)
      juce::SpinLock::ScopedLockType profileLockGuard(profileLock);

      bool hasCh1 = false;
      for (int i = 0; i < profilesTree.getNumChildren(); ++i) {
        juce::ValueTree node = profilesTree.getChild(i);
        int channel = node.getProperty("channel", 0);
        if (channel == 1)
          hasCh1 = true;

        int mode = node.getProperty("mode", 0);
        uint32_t size = static_cast<uint32_t>(
            static_cast<int>(node.getProperty("size", 0)));
        uint32_t blockCount = static_cast<uint32_t>(
            static_cast<int>(node.getProperty("blockCount", 0)));
        juce::String base64Data = node.getProperty("data", "");

        if (size > 0 && base64Data.isNotEmpty() && mode >= 1 && mode <= 4) {
          juce::MemoryBlock mb;
          if (mb.fromBase64Encoding(base64Data) &&
              mb.getSize() >= size * sizeof(float)) {
            const float* floatArray =
                reinterpret_cast<const float*>(mb.getData());

            // Narrow (or not yet created) engines cannot take this channel —
            // keep it pending for a later, wider rebuild instead of dropping
            // it.
            if (!loadProfileResampledChecked(
                    engineGroup.get(), channel, floatArray, size, blockCount,
                    static_cast<SpecbleachProfileMode>(mode))) {
              PendingProfile pp;
              pp.channel = channel;
              pp.mode = mode;
              pp.size = size;
              pp.blockCount = blockCount;
              pp.data.assign(floatArray, floatArray + size);
              pendingProfiles.push_back(std::move(pp));
            }
          }
        }
      }

      // If legacy profile had no Channel 1 profile, copy Channel 0 to
      // Channel 1 (stereo groups only)
      if (!hasCh1 &&
          specbleach_stereo_get_channel_count(engineGroup.get()) > 1) {
        for (const auto& pp : pendingProfiles) {
          if (pp.channel == 0) {
            loadProfileResampledChecked(
                engineGroup.get(), 1, pp.data.data(), pp.size, pp.blockCount,
                static_cast<SpecbleachProfileMode>(pp.mode));
          }
        }
      }

      // Fill any remaining per-channel/per-mode gaps
      specbleach_stereo_sync_profiles(engineGroup.get());
    }
  }
}

juce::AudioProcessorEditor* NoiseRepellentAudioProcessor::createEditor() {
  return new NoiseRepellentAudioProcessorEditor(*this);
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
  return new NoiseRepellentAudioProcessor();
}
