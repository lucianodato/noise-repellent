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

NoiseRepellentAudioProcessor::NoiseRepellentAudioProcessor()
    : AudioProcessor(
          BusesProperties()
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      parameters(*this, nullptr, "PARAMETERS", createParameterLayout()),
      dryWetMixer(16384) {
  bypassParameter = dynamic_cast<juce::AudioParameterBool*>(
      parameters.getParameter("bypass"));
  startTimerHz(60); // deferred latency reporting to the message thread
  DBG("libspecbleach " << specbleach_get_version_string());
}

NoiseRepellentAudioProcessor::~NoiseRepellentAudioProcessor() {
  stopTimer();
  releaseResources();
}

juce::AudioProcessorParameter*
NoiseRepellentAudioProcessor::getBypassParameter() const {
  return bypassParameter;
}

juce::AudioProcessorValueTreeState::ParameterLayout
NoiseRepellentAudioProcessor::createParameterLayout() {
  std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

  params.push_back(std::make_unique<juce::AudioParameterChoice>(
      "algorithm_mode", "Algorithm",
      juce::StringArray{"1D Spectral", "2D NLM Patch HQ"}, 0));

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
    m.front() = (dX.front() > 1e-5f) ? dY.front() / dX.front() : 0.0f;
    m.back() = (dX.back() > 1e-5f) ? dY.back() / dX.back() : 0.0f;

    for (size_t i = 1; i < numNodes - 1; ++i) {
      float secant1 = (dX[i - 1] > 1e-5f) ? dY[i - 1] / dX[i - 1] : 0.0f;
      float secant2 = (dX[i] > 1e-5f) ? dY[i] / dX[i] : 0.0f;
      if (secant1 * secant2 <= 0.0f) {
        m[i] = 0.0f;
      } else {
        m[i] = 0.5f * (secant1 + secant2);
      }
    }
  } else {
    float secant = (dX.front() > 1e-5f) ? dY.front() / dX.front() : 0.0f;
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
        float t = (dx > 1e-9f) ? (binNormX - sortedNodes[i].normX) / dx : 0.0f;
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

// Resamples a persisted noise profile to whatever spectrum size the target
// group expects, then loads it into one of its channel engines.
static bool loadProfileGroupSafe(specbleach_stereo* group, uint32_t channel,
                                 const float* data, uint32_t size,
                                 uint32_t blockCount, int mode) {
  if (group == nullptr || data == nullptr || size == 0)
    return false;

  const uint32_t targetSize = specbleach_stereo_get_noise_profile_size(group);
  if (targetSize == 0)
    return false;

  if (size == targetSize) {
    return specbleach_stereo_load_noise_profile_for_channel(
        group, channel, data, size, blockCount, mode);
  } else if (targetSize == 1) {
    float resampledSample = data[0];
    return specbleach_stereo_load_noise_profile_for_channel(
        group, channel, &resampledSample, 1, blockCount, mode);
  } else {
    std::vector<float> resampled(targetSize);
    float maxSrcIdx = static_cast<float>(size - 1);
    float maxTargetIdx = static_cast<float>(targetSize - 1);
    for (uint32_t i = 0; i < targetSize; ++i) {
      float normPos = static_cast<float>(i) / maxTargetIdx;
      float exactIdx = normPos * maxSrcIdx;
      uint32_t idx0 = std::clamp(static_cast<uint32_t>(exactIdx), 0u, size - 1);
      uint32_t idx1 = std::min(idx0 + 1, size - 1);
      float frac = exactIdx - static_cast<float>(idx0);
      resampled[i] = (1.0f - frac) * data[idx0] + frac * data[idx1];
    }
    return specbleach_stereo_load_noise_profile_for_channel(
        group, channel, resampled.data(), targetSize, blockCount, mode);
  }
}

specbleach_stereo* NoiseRepellentAudioProcessor::activeGroupFor(
    int algoMode) const {
  return (algoMode == 1) ? nlmGroup.get() : spectralGroup.get();
}

void NoiseRepellentAudioProcessor::ensureEnginesInitialized(double sampleRate) {
  bool needNewEngines = (spectralGroup == nullptr || nlmGroup == nullptr ||
                         std::abs(currentSampleRate - sampleRate) > 0.001);

  if (!needNewEngines)
    return;

  struct SavedProfileData {
    int channel; // 0 = Left, 1 = Right
    int mode;
    uint32_t size;
    uint32_t blockCount;
    std::vector<float> data;
  };

  std::vector<SavedProfileData> profilesSpectral;
  std::vector<SavedProfileData> profilesNLM;

  auto backupGroup = [](specbleach_stereo* group,
                        std::vector<SavedProfileData>& out) {
    if (group == nullptr)
      return;
    const uint32_t channels = specbleach_stereo_get_channel_count(group);
    const uint32_t sz = specbleach_stereo_get_noise_profile_size(group);
    for (uint32_t channel = 0; channel < channels; ++channel) {
      for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
           mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
        if (!specbleach_stereo_profile_available_for_channel(group, channel,
                                                             mode))
          continue;
        float* profile = specbleach_stereo_get_noise_profile_for_channel(
            group, channel, mode);
        const uint32_t blockCount =
            specbleach_stereo_get_profile_block_count_for_channel(
                group, channel, mode);
        if (profile != nullptr && sz > 0) {
          out.push_back({static_cast<int>(channel), mode, sz, blockCount,
                         std::vector<float>(profile, profile + sz)});
        }
      }
    }
  };

  backupGroup(spectralGroup.get(), profilesSpectral);
  backupGroup(nlmGroup.get(), profilesNLM);

  // Free existing groups before allocating replacements
  spectralGroup.reset();
  nlmGroup.reset();

  uint32_t channels = preparedNumChannels;
  if (channels < 2) {
    channels = static_cast<uint32_t>(
        std::max({getTotalNumInputChannels(), getTotalNumOutputChannels(), 2}));
  }

  const uint32_t sampleRateUint = static_cast<uint32_t>(sampleRate);
  spectralGroup = specbleach::make_stereo_group(
      sampleRateUint, 50.0f, channels, SPECBLEACH_STEREO_ENGINE_SPECTRAL);
  nlmGroup = specbleach::make_stereo_group(sampleRateUint, 50.0f, channels,
                                           SPECBLEACH_STEREO_ENGINE_NLM_2D);

  currentSampleRate = sampleRate;
  preparedNumChannels = channels;

  // Restore backed-up profiles into the new groups
  for (const auto& item : profilesSpectral) {
    loadProfileGroupSafe(spectralGroup.get(),
                         static_cast<uint32_t>(item.channel), item.data.data(),
                         item.size, item.blockCount, item.mode);
  }

  for (const auto& item : profilesNLM) {
    loadProfileGroupSafe(nlmGroup.get(), static_cast<uint32_t>(item.channel),
                         item.data.data(), item.size, item.blockCount,
                         item.mode);
  }

  // Restore state pending profiles if loaded before prepareToPlay
  bool hasCh1 = false;
  for (const auto& item : pendingProfiles) {
    if (item.channel == 1)
      hasCh1 = true;

    loadProfileGroupSafe(spectralGroup.get(),
                         static_cast<uint32_t>(item.channel), item.data.data(),
                         item.size, item.blockCount, item.mode);
    loadProfileGroupSafe(nlmGroup.get(), static_cast<uint32_t>(item.channel),
                         item.data.data(), item.size, item.blockCount,
                         item.mode);
  }

  // Legacy states may lack channel 1 profiles — fall back from channel 0
  if (!hasCh1) {
    for (const auto& item : pendingProfiles) {
      if (item.channel == 0) {
        loadProfileGroupSafe(spectralGroup.get(), 1u, item.data.data(),
                             item.size, item.blockCount, item.mode);
        loadProfileGroupSafe(nlmGroup.get(), 1u, item.data.data(), item.size,
                             item.blockCount, item.mode);
      }
    }
  }
  pendingProfiles.clear();

  // Fill any remaining per-channel/per-mode gaps within each group
  specbleach_stereo_sync_profiles(spectralGroup.get());
  specbleach_stereo_sync_profiles(nlmGroup.get());
}

void NoiseRepellentAudioProcessor::prepareToPlay(double sampleRate,
                                                 int samplesPerBlock) {
  ensureEnginesInitialized(sampleRate);

  currentAlgoMode = static_cast<int>(
      parameters.getRawParameterValue("algorithm_mode")->load());

  uint32_t nativeLatency = 0;
  if (auto* activeGroup = activeGroupFor(currentAlgoMode))
    nativeLatency = specbleach_stereo_get_latency(activeGroup);
  setLatencySamples(static_cast<int>(nativeLatency));

  const int bufferCapacity = std::max(samplesPerBlock, 16384);
  preparedBlockSize = static_cast<uint32_t>(bufferCapacity);
  preparedNumChannels = static_cast<uint32_t>(
      std::max({getTotalNumInputChannels(), getTotalNumOutputChannels(), 2}));

  juce::dsp::ProcessSpec spec;
  spec.sampleRate = sampleRate;
  spec.maximumBlockSize = preparedBlockSize;
  spec.numChannels = preparedNumChannels;

  dryWetMixer.prepare(spec);
  dryWetMixer.setMixingRule(juce::dsp::DryWetMixingRule::linear);
  dryWetMixer.setWetLatency(static_cast<float>(nativeLatency));

  currentLatency = static_cast<int>(nativeLatency);
  visualizerDelayBuffer.assign(
      std::max<size_t>(32768, static_cast<size_t>(nativeLatency) + 8192), 0.0f);
  visualizerDelayWritePos = 0;

  // Pre-allocate persistent buffers to prevent audio-thread allocations
  dryInputL.resize(static_cast<size_t>(bufferCapacity), 0.0f);
  wetScratchA.setSize(static_cast<int>(preparedNumChannels), bufferCapacity,
                      false, false, true);
  wetScratchB.setSize(static_cast<int>(preparedNumChannels), bufferCapacity,
                      false, false, true);
  jassert(dryInputL.size() >= static_cast<size_t>(samplesPerBlock));
  jassert(wetScratchA.getNumSamples() >= samplesPerBlock);

  switchPhase = SwitchPhase::Steady;
  stageSamplesRemaining = 0;
  warmSamplesTotal = 0;
  latencyAnnounceCountdown = -1;
  muteGain = 1.0f;
  uiSwitchProgress.store(1.0f, std::memory_order_relaxed);

  currentLatency = static_cast<int>(nativeLatency);
  dryWetMixer.setWetLatency(static_cast<float>(nativeLatency));

  // Reset FFT accumulation
  fftAccumInput.fill(0.0f);
  fftAccumOutput.fill(0.0f);
  fftAccumTransient.fill(0.0f);
  fftAccumCount = 0;
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

void NoiseRepellentAudioProcessor::processBlock(
    juce::AudioBuffer<float>& buffer, juce::MidiBuffer&) {
  juce::ScopedNoDenormals noDenormals;
  const int numSamples = buffer.getNumSamples();
  const int numChannels = buffer.getNumChannels();

  if (numSamples == 0 || numChannels == 0)
    return;

  const bool isBypassed =
      parameters.getRawParameterValue("bypass")->load() > 0.5f;
  const int algoMode = static_cast<int>(
      parameters.getRawParameterValue("algorithm_mode")->load());
  const bool transientProtectionEnable =
      parameters.getRawParameterValue("transient_protection_enable")->load() >
      0.5f;
  const bool learnNoise =
      parameters.getRawParameterValue("learn_noise")->load() > 0.5f;
  const bool adaptiveNoise =
      parameters.getRawParameterValue("adaptive_noise")->load() > 0.5f;
  const int adaptiveMethod = static_cast<int>(
      parameters.getRawParameterValue("adaptive_method")->load());
  const float aggressiveness =
      parameters.getRawParameterValue("aggressiveness")->load();

  const bool linkReduction =
      parameters.getRawParameterValue("link_reduction")->load() > 0.5f;
  const float masterRed =
      parameters.getRawParameterValue("reduction_amount")->load();
  const float tonalRed =
      linkReduction
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
      linkThresholdOffset
          ? profileOffset
          : parameters.getRawParameterValue("tonal_noise_profile_offset")
                ->load();
  const bool curveEnabled =
      parameters.getRawParameterValue("reduction_curve_enabled")->load() > 0.5f;

  if (curveEnabled) {
    if (curveNodesDirty.exchange(false)) {
      juce::SpinLock::ScopedTryLockType tryLock(curveLock);
      if (tryLock.isLocked()) {
        uint32_t profileSize = 0;
        if (auto* activeGrp = activeGroupFor(algoMode))
          profileSize = specbleach_stereo_get_noise_profile_size(activeGrp);
        if (profileSize > 0)
          interpolateCurve(profileSize);
      } else {
        curveNodesDirty = true;
      }
    }
  }

  // Migrate profiles between engine groups when manual learning was just
  // turned off. Inline on the audio thread like the original plugin: the
  // migration is a bounded, allocation-free copy of small profile buffers.
  if (wasLearning && !learnNoise && spectralGroup != nullptr &&
      nlmGroup != nullptr) {
    auto* sourceGroup = activeGroupFor(currentAlgoMode);
    auto* otherGroup =
        (currentAlgoMode == 1) ? spectralGroup.get() : nlmGroup.get();
    if (sourceGroup != nullptr && otherGroup != nullptr &&
        sourceGroup != otherGroup) {
      specbleach_stereo_migrate_profiles_from(otherGroup, sourceGroup);
      specbleach_stereo_sync_profiles(otherGroup);
    }
  }
  wasLearning = learnNoise;

  // Compute linear parameters for libspecbleach
  const float reductionGain = juce::Decibels::decibelsToGain(-masterRed);
  const float tonalReductionGain = juce::Decibels::decibelsToGain(-tonalRed);
  const float smoothingNorm = smoothing / 100.0f;
  const float whiteningNorm = whitening / 100.0f;
  const float profileScale = std::pow(10.0f, profileOffset / 10.0f);
  const float tonalProfileScale = std::pow(10.0f, tonalProfileOffsetDb / 10.0f);

  // Prepare parameters for DSP engines
  SpecbleachDenoiserParameters p{};
  p.learn_noise = learnNoise ? SPECBLEACH_LEARN_ALL : SPECBLEACH_LEARN_OFF;
  p.residual_listen = residualListen;
  p.reduction_gain = reductionGain;
  p.smoothing_factor = smoothingNorm;
  p.whitening_factor = whiteningNorm;
  p.adaptive_noise = adaptiveNoise;
  p.noise_estimation_method =
      static_cast<SpecbleachNoiseEstimationMethod>(adaptiveMethod);
  p.masking_depth = masking;
  p.suppression_strength = 1.0f;
  p.aggressiveness = aggressiveness;
  p.tonal_reduction_gain = tonalReductionGain;
  p.hpss_enable = transientProtectionEnable;
  p.noise_profile_scale = profileScale;
  p.reduction_curve_bias =
      curveEnabled ? interpolatedCurveBias.data() : nullptr;
  p.reduction_curve_enabled = curveEnabled;
  p.reduction_curve_size =
      curveEnabled ? static_cast<uint32_t>(interpolatedCurveBias.size()) : 0;
  p.tonal_noise_profile_scale = tonalProfileScale;

  Specbleach2DDenoiserParameters p2{};
  p2.learn_noise = learnNoise ? SPECBLEACH_LEARN_ALL : SPECBLEACH_LEARN_OFF;
  p2.residual_listen = residualListen;
  p2.reduction_gain = reductionGain;
  p2.smoothing_factor = smoothingNorm;
  p2.whitening_factor = whiteningNorm;
  p2.adaptive_noise = adaptiveNoise;
  p2.noise_estimation_method =
      static_cast<SpecbleachNoiseEstimationMethod>(adaptiveMethod);
  p2.nlm_masking_protection = masking;
  p2.suppression_strength = 1.0f;
  p2.aggressiveness = aggressiveness;
  p2.tonal_reduction_gain = tonalReductionGain;
  p2.hpss_enable = transientProtectionEnable;
  p2.noise_profile_scale = profileScale;
  p2.reduction_curve_bias =
      curveEnabled ? interpolatedCurveBias.data() : nullptr;
  p2.reduction_curve_enabled = curveEnabled;
  p2.reduction_curve_size =
      curveEnabled ? static_cast<uint32_t>(interpolatedCurveBias.size()) : 0;
  p2.tonal_noise_profile_scale = tonalProfileScale;

  // ── Engine switch: fade out -> silent warm-up -> fade in ─────────────
  const int edgeFadeSamples =
      static_cast<int>(currentSampleRate * (kEdgeFadeMs / 1000.0));

  if (algoMode != currentAlgoMode && switchPhase == SwitchPhase::Steady &&
      spectralGroup != nullptr && nlmGroup != nullptr) {
    auto* sourceGroup = activeGroupFor(currentAlgoMode);
    auto* targetGroup = activeGroupFor(algoMode);
    if (sourceGroup != nullptr && targetGroup != nullptr &&
        sourceGroup != targetGroup) {
      specbleach_stereo_migrate_profiles_from(targetGroup, sourceGroup);
      specbleach_stereo_sync_profiles(targetGroup);

      switchFromMode = currentAlgoMode;
      currentAlgoMode = algoMode;
      switchPhase = SwitchPhase::FadeOut;
      stageSamplesRemaining = edgeFadeSamples;
      warmSamplesTotal =
          static_cast<int>(currentSampleRate * (kWarmupMs / 1000.0));
      latencyAnnounceCountdown =
          static_cast<int>(currentSampleRate * (kLatencySettleMs / 1000.0));
      // NOTE: reported latency moves at the FadeOut->WarmSilent boundary,
      // once the output is fully silent. Announcing a latency DROP while
      // audio is still fading makes hosts truncate in-flight compensation
      // and glitch the tail (the 2D->1D artifact). The handoff itself is
      // staged to the message thread (see pendingLatencyReport): hosts like
      // Reaper suspend/resume the plugin on synchronous latency changes.
    }
  }

  const bool switchingNow = (switchPhase != SwitchPhase::Steady);

  // Save dry input copy for FFT visualization before processing
  const size_t copySamples =
      std::min(static_cast<size_t>(numSamples), dryInputL.size());
  if (numChannels >= 1 && copySamples > 0)
    std::copy_n(buffer.getReadPointer(0), copySamples, dryInputL.begin());

  // Push dry samples into JUCE DryWetMixer before in-place DSP processing
  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);

  if (!switchingNow) {
    // Steady state: exactly one family, in place, native latency
    auto* activeGroup = activeGroupFor(currentAlgoMode);
    if (activeGroup != nullptr) {
      const uint32_t groupChannels = std::min<uint32_t>(
          specbleach_stereo_get_channel_count(activeGroup), 2u);
      const float* inPtrs[2] = {nullptr, nullptr};
      float* sinkA = wetScratchA.getNumChannels() > 1
                         ? wetScratchA.getWritePointer(1)
                         : nullptr;
      for (uint32_t ch = 0; ch < groupChannels; ++ch)
        inPtrs[ch] = buffer.getReadPointer(
            ch < static_cast<uint32_t>(numChannels) ? static_cast<int>(ch) : 0);

      float* outPtrs[2] = {nullptr, nullptr};
      for (uint32_t ch = 0; ch < groupChannels; ++ch) {
        outPtrs[ch] =
            (ch < static_cast<uint32_t>(numChannels))
                ? buffer.getWritePointer(static_cast<int>(ch))
                : ((sinkA != nullptr) ? sinkA : buffer.getWritePointer(0));
      }

      if (currentAlgoMode == 0 && activeGroup == spectralGroup.get()) {
        specbleach_stereo_load_parameters_1d(spectralGroup.get(), &p,
                                             sizeof(p));
      } else if (currentAlgoMode == 1 && activeGroup == nlmGroup.get()) {
        specbleach_stereo_load_parameters_2d(nlmGroup.get(), &p2, sizeof(p2));
      }
      specbleach_stereo_process(activeGroup, static_cast<uint32_t>(numSamples),
                                inPtrs, outPtrs);
    }
  } else {
    // Switching stages: FadeOut keeps the CURRENT family audible while its
    // gain ramps to zero; WarmSilent renders only the target (discarded) so
    // it is warm when audio returns; FadeIn ramps the new family up.
    int renderMode;
    if (switchPhase == SwitchPhase::FadeOut)
      renderMode = switchFromMode;
    else
      renderMode = currentAlgoMode;

    auto* renderGroup = activeGroupFor(renderMode);
    if (renderGroup != nullptr) {
      const uint32_t groupChannels = std::min<uint32_t>(
          specbleach_stereo_get_channel_count(renderGroup), 2u);
      const float* inPtrs[2] = {nullptr, nullptr};
      for (uint32_t ch = 0; ch < groupChannels; ++ch)
        inPtrs[ch] = buffer.getReadPointer(
            ch < static_cast<uint32_t>(numChannels) ? static_cast<int>(ch) : 0);

      float* sinkA = wetScratchA.getNumChannels() > 1
                         ? wetScratchA.getWritePointer(1)
                         : nullptr;
      float* outPtrs[2] = {nullptr, nullptr};
      for (uint32_t ch = 0; ch < groupChannels; ++ch) {
        outPtrs[ch] =
            (ch < static_cast<uint32_t>(numChannels))
                ? buffer.getWritePointer(static_cast<int>(ch))
                : ((sinkA != nullptr) ? sinkA : buffer.getWritePointer(0));
      }

      if (renderMode == 0 && renderGroup == spectralGroup.get()) {
        specbleach_stereo_load_parameters_1d(spectralGroup.get(), &p,
                                             sizeof(p));
      } else if (renderMode == 1 && renderGroup == nlmGroup.get()) {
        specbleach_stereo_load_parameters_2d(nlmGroup.get(), &p2, sizeof(p2));
      }
      specbleach_stereo_process(renderGroup, static_cast<uint32_t>(numSamples),
                                inPtrs, outPtrs);

      if (switchPhase == SwitchPhase::WarmSilent)
        buffer.clear(); // nothing audible during the silent window
    }

    // Latency handoff: only after the old tail drained through the host;
    // delivered asynchronously by the timer
    if (switchPhase == SwitchPhase::WarmSilent &&
        latencyAnnounceCountdown > 0) {
      latencyAnnounceCountdown -= numSamples;
      if (latencyAnnounceCountdown <= 0) {
        latencyAnnounceCountdown = -1;
        if (auto* targetGroup = activeGroupFor(currentAlgoMode)) {
          pendingLatencyReport.store(
              static_cast<int>(specbleach_stereo_get_latency(targetGroup)),
              std::memory_order_release);
        }
      }
    }

    // Gain envelope: linear toward the stage's destination across the edge
    const float gainTarget =
        (switchPhase == SwitchPhase::WarmSilent) ? 0.0f : 1.0f;
    const float step = 1.0f / static_cast<float>(std::max(1, edgeFadeSamples));
    const int nsamp = numSamples;
    if (muteGain != gainTarget || gainTarget == 0.0f) {
      for (int chI = 0; chI < buffer.getNumChannels(); ++chI) {
        float* d = buffer.getWritePointer(chI);
        float g = muteGain;
        for (int smp = 0; smp < nsamp; ++smp) {
          if (gainTarget < g)
            g = std::max(0.0f, g - step);
          else
            g = std::min(1.0f, g + step);
          d[smp] *= g;
        }
      }
      muteGain = gainTarget;
    }

    // Stage advancement
    stageSamplesRemaining -= numSamples;
    if (stageSamplesRemaining <= 0) {
      switch (switchPhase) {
        case SwitchPhase::FadeOut:
          switchPhase = SwitchPhase::WarmSilent;
          stageSamplesRemaining = warmSamplesTotal;
          break;
        case SwitchPhase::WarmSilent:
          switchPhase = SwitchPhase::FadeIn;
          stageSamplesRemaining = edgeFadeSamples;
          break;
        case SwitchPhase::FadeIn:
          switchPhase = SwitchPhase::Steady;
          break;
        default:
          break;
      }
    }
  }

  // UI progress across the three switch stages
  {
    float uiP = 1.0f;
    if (switchPhase != SwitchPhase::Steady) {
      const float total =
          static_cast<float>(2 * edgeFadeSamples + warmSamplesTotal);
      float elapsed = 0.0f;
      if (switchPhase == SwitchPhase::FadeOut) {
        elapsed = static_cast<float>(edgeFadeSamples - stageSamplesRemaining);
      } else if (switchPhase == SwitchPhase::WarmSilent) {
        elapsed = static_cast<float>(edgeFadeSamples) +
                  static_cast<float>(warmSamplesTotal - stageSamplesRemaining);
      } else {
        elapsed = static_cast<float>(2 * edgeFadeSamples + warmSamplesTotal -
                                     stageSamplesRemaining);
      }
      uiP = total > 0.0f ? elapsed / total : 1.0f;
    }
    uiSwitchProgress.store(juce::jlimit(0.0f, 1.0f, uiP),
                           std::memory_order_relaxed);
  }

  // Query transient protection status and intensity (aggregated max across
  // channels; report the TARGET family while a fade is running)
  float reportedTransientIntensity = 0.0f;
  {
    auto* reportGroup = activeGroupFor(currentAlgoMode);
    if (reportGroup != nullptr)
      reportedTransientIntensity =
          specbleach_stereo_get_transient_intensity(reportGroup);
  }
  transientActivity.store(reportedTransientIntensity,
                          std::memory_order_relaxed);
  transientProtectionActive.store(transientProtectionEnable,
                                  std::memory_order_relaxed);

  // Apply Soft Crossfade Bypass using JUCE DryWetMixer (with dry latency
  // compensation)
  dryWetMixer.setWetMixProportion(isBypassed ? 0.0f : 1.0f);
  dryWetMixer.mixWetSamples(audioBlock);

  // Skip FFT analysis and FIFO writes during offline rendering or if the GUI is
  // closed
  if (isNonRealtime() || getActiveEditor() == nullptr) {
    fftAccumCount = 0;
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

        bool isAdaptive = adaptiveNoise;
        bool isLearning = learnNoise;
        bool profileHasAnyMode = false;

        juce::SpinLock::ScopedTryLockType profileTryLock(profileLock);
        if (profileTryLock.isLocked()) {
          auto* vizGroup = activeGroupFor(algoMode);
          if (vizGroup != nullptr) {
            for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
                 mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
              if (specbleach_stereo_profile_available_for_channel(vizGroup, 0u,
                                                                  mode)) {
                profileHasAnyMode = true;
                break;
              }
            }

            if (isLearning) {
              actualNoiseProfile =
                  specbleach_stereo_get_noise_profile_for_channel(vizGroup, 0u,
                                                                  1);
              profileSize = specbleach_stereo_get_noise_profile_size(vizGroup);
              profileAvailable =
                  (actualNoiseProfile != nullptr && profileSize > 0);
            } else if (isAdaptive || profileHasAnyMode) {
              actualNoiseProfile =
                  specbleach_stereo_get_active_noise_profile_for_channel(
                      vizGroup, 0u);
              profileSize = specbleach_stereo_get_noise_profile_size(vizGroup);
              if (!actualNoiseProfile && profileHasAnyMode) {
                for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
                     mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
                  if (specbleach_stereo_profile_available_for_channel(
                          vizGroup, 0u, mode)) {
                    actualNoiseProfile =
                        specbleach_stereo_get_noise_profile_for_channel(
                            vizGroup, 0u, mode);
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
              (parameters.getRawParameterValue("link_reduction")->load() >
               0.5f);
          frame.isOffsetLinked =
              (parameters.getRawParameterValue("link_threshold_offset")
                   ->load() > 0.5f);
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
          frame.isTransientProtectionActive = transientProtectionEnable;
          frame.tonalPeaksHz.clear();

          if (profileAvailable && actualNoiseProfile) {
            // profileSize represents the unique spectrum from DC (0 Hz) to
            // Nyquist (Fs/2)
            size_t realProfileBins = profileSize;

            const float maxProfileIdx = static_cast<float>(realProfileBins - 1);
            const float maxFftIdx = static_cast<float>(kFftBins - 1);
            // FFTW unnormalized power scaling offset: 20 * log10(N/2)
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

              float val0 = actualNoiseProfile[p0];
              float val1 = actualNoiseProfile[p1];
              float interpVal = (1.0f - frac) * val0 + frac * val1;

              float rawDb = 10.0f * std::log10(std::max(interpVal, 1e-12f));
              frame.noiseFloorDB[i] = rawDb - dbOffset;
            }

            // Query libspecbleach directly for reported tonal peak frequencies
            // computed on this real noise profile
            std::array<float, 32> peakBuf{};
            uint32_t numPeaks = 0;
            if (auto* vizGroup = activeGroupFor(algoMode)) {
              numPeaks = specbleach_stereo_get_tonal_peaks_for_channel(
                  vizGroup, 0u, peakBuf.data(),
                  static_cast<uint32_t>(peakBuf.size()));
            }

            for (uint32_t i = 0; i < numPeaks; ++i) {
              frame.tonalPeaksHz.push_back(peakBuf[i]);
            }
          } else {
            frame.noiseFloorDB.fill(-120.0f);
          }
        }

        spectralFifo.finishedWrite(1);
      }

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
  specbleach_stereo_reset_profiles(spectralGroup.get());
  specbleach_stereo_reset_profiles(nlmGroup.get());
  pendingProfiles.clear();
}

bool NoiseRepellentAudioProcessor::hasNoiseProfile() const {
  for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
       mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
    if (spectralGroup != nullptr &&
        specbleach_stereo_profile_available_for_channel(spectralGroup.get(), 0u,
                                                        mode))
      return true;
    if (nlmGroup != nullptr && specbleach_stereo_profile_available_for_channel(
                                   nlmGroup.get(), 0u, mode))
      return true;
  }
  return false;
}

void NoiseRepellentAudioProcessor::timerCallback() {
  const int staged =
      pendingLatencyReport.exchange(-1, std::memory_order_acq_rel);
  if (staged >= 0) {
    setLatencySamples(staged);
    dryWetMixer.setWetLatency(static_cast<float>(staged));
  }
}

juce::AudioProcessorEditor* NoiseRepellentAudioProcessor::createEditor() {
  return new NoiseRepellentAudioProcessorEditor(*this);
}

void NoiseRepellentAudioProcessor::getStateInformation(
    juce::MemoryBlock& destData) {
  juce::SpinLock::ScopedLockType lock(profileLock);

  // Fill per-channel/per-mode gaps before serialization
  specbleach_stereo_sync_profiles(spectralGroup.get());
  specbleach_stereo_sync_profiles(nlmGroup.get());

  auto state = parameters.copyState();

  juce::ValueTree profilesTree("LEARNED_PROFILES");

  const uint32_t channelCount =
      (spectralGroup != nullptr)
          ? specbleach_stereo_get_channel_count(spectralGroup.get())
          : ((nlmGroup != nullptr)
                 ? specbleach_stereo_get_channel_count(nlmGroup.get())
                 : 0);

  for (uint32_t channel = 0; channel < channelCount; ++channel) {
    for (int mode = SPECBLEACH_PROFILE_MODE_FIRST;
         mode <= SPECBLEACH_PROFILE_MODE_LAST; ++mode) {
      const float* profile = nullptr;
      uint32_t profileSize = 0;
      uint32_t blockCount = 0;

      auto collectFrom = [&](specbleach_stereo* group) {
        if (group == nullptr || profile != nullptr)
          return;
        if (specbleach_stereo_profile_available_for_channel(group, channel,
                                                            mode)) {
          profile = specbleach_stereo_get_noise_profile_for_channel(
              group, channel, mode);
          profileSize = specbleach_stereo_get_noise_profile_size(group);
          blockCount = specbleach_stereo_get_profile_block_count_for_channel(
              group, channel, mode);
        }
      };

      collectFrom(spectralGroup.get());
      collectFrom(nlmGroup.get());

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

            PendingProfile pp;
            pp.channel = channel;
            pp.mode = mode;
            pp.size = size;
            pp.blockCount = blockCount;
            pp.data.assign(floatArray, floatArray + size);
            pendingProfiles.push_back(pp);

            const uint32_t channelUint = static_cast<uint32_t>(channel);
            loadProfileGroupSafe(spectralGroup.get(), channelUint, floatArray,
                                 size, blockCount, mode);
            loadProfileGroupSafe(nlmGroup.get(), channelUint, floatArray, size,
                                 blockCount, mode);
          }
        }
      }

      // If legacy profile had no Channel 1 profile, copy Channel 0 to Channel 1
      // in both groups
      if (!hasCh1) {
        for (const auto& pp : pendingProfiles) {
          if (pp.channel == 0) {
            loadProfileGroupSafe(spectralGroup.get(), 1u, pp.data.data(),
                                 pp.size, pp.blockCount, pp.mode);
            loadProfileGroupSafe(nlmGroup.get(), 1u, pp.data.data(), pp.size,
                                 pp.blockCount, pp.mode);
          }
        }
      }

      // Fill any remaining per-channel/per-mode gaps within each group
      specbleach_stereo_sync_profiles(spectralGroup.get());
      specbleach_stereo_sync_profiles(nlmGroup.get());

      if (spectralGroup != nullptr && nlmGroup != nullptr) {
        pendingProfiles.clear();
      }
    }
  }
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
  return new NoiseRepellentAudioProcessor();
}
