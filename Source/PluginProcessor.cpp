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
}

NoiseRepellentAudioProcessor::~NoiseRepellentAudioProcessor() {
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
      juce::StringArray{"1D Spectral", "2D NLM Patch HQ"}, 1));

  params.push_back(std::make_unique<juce::AudioParameterChoice>(
      "hpss_quality", "HPSS Quality",
      juce::StringArray{"Off", "Low (Fast)", "Mid (Balanced)", "High (Ultra)"},
      0));

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
      juce::NormalisableRange<float>(-1.0f, 1.0f, 0.1f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "link_reduction", "Link Reduction", true,
      juce::AudioParameterBoolAttributes().withAutomatable(false).withMeta(
          true)));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "reduction_amount", "Master Reduction",
      juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 15.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "tonal_reduction", "Tonal Reduction",
      juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 15.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "smoothing_factor", "Smoothing",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "masking_depth", "Masking Protect",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 100.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "whitening_factor", "Whitening",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "suppression_strength", "Suppression",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 50.0f));

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
              -sortedNodes.front().biasDB);
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
      interpolatedCurveBias[k] = -sortedNodes.front().biasDB;
      continue;
    }
    if (binNormX >= sortedNodes.back().normX) {
      interpolatedCurveBias[k] = -sortedNodes.back().biasDB;
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
        interpolatedCurveBias[k] = -val;
        break;
      }
    }
  }
}

template <typename SizeQueryFn, typename LoadProfileFn>
static bool loadProfileSafeHelper(SpectralBleachHandle h, const float* data,
                                  uint32_t size, uint32_t blockCount, int mode,
                                  SizeQueryFn sizeQueryFn,
                                  LoadProfileFn loadProfileFn) {
  if (!h || !data || size == 0)
    return false;
  uint32_t targetSize = sizeQueryFn(h);
  if (targetSize == 0)
    return false;

  if (size == targetSize) {
    return loadProfileFn(h, data, size, blockCount, mode);
  } else if (targetSize == 1) {
    float resampledSample = data[0];
    return loadProfileFn(h, &resampledSample, 1, blockCount, mode);
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
    return loadProfileFn(h, resampled.data(), targetSize, blockCount, mode);
  }
}

static bool loadProfile1DSafe(SpectralBleachHandle h, const float* data,
                              uint32_t size, uint32_t blockCount, int mode) {
  return loadProfileSafeHelper(h, data, size, blockCount, mode,
                               specbleach_get_noise_profile_size,
                               specbleach_load_noise_profile_for_mode);
}

static bool loadProfile2DSafe(SpectralBleachHandle h, const float* data,
                              uint32_t size, uint32_t blockCount, int mode) {
  return loadProfileSafeHelper(h, data, size, blockCount, mode,
                               specbleach_2d_get_noise_profile_size,
                               specbleach_2d_load_noise_profile_for_mode);
}

void NoiseRepellentAudioProcessor::ensureEnginesInitialized(double sampleRate) {
  bool needNewEngines = (specbleach1D_L == nullptr ||
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

  std::vector<SavedProfileData> profiles1D;
  std::vector<SavedProfileData> profiles2D;

  auto backup1D = [&](SpectralBleachHandle h, int channel) {
    if (!h)
      return;
    for (int mode = 1; mode <= 4; ++mode) {
      if (specbleach_noise_profile_available_for_mode(h, mode)) {
        float* p = specbleach_get_noise_profile_for_mode(h, mode);
        uint32_t sz = specbleach_get_noise_profile_size(h);
        uint32_t bc =
            specbleach_get_noise_profile_block_count_for_mode(h, mode);
        if (p && sz > 0) {
          profiles1D.push_back(
              {channel, mode, sz, bc, std::vector<float>(p, p + sz)});
        }
      }
    }
  };
  backup1D(specbleach1D_L, 0);
  backup1D(specbleach1D_R, 1);

  auto backup2D = [&](auto h, int channel) {
    if (!h)
      return;
    for (int mode = 1; mode <= 4; ++mode) {
      if (specbleach_2d_noise_profile_available_for_mode(h, mode)) {
        float* p = specbleach_2d_get_noise_profile_for_mode(h, mode);
        uint32_t sz = specbleach_2d_get_noise_profile_size(h);
        uint32_t bc =
            specbleach_2d_get_noise_profile_block_count_for_mode(h, mode);
        if (p && sz > 0) {
          profiles2D.push_back(
              {channel, mode, sz, bc, std::vector<float>(p, p + sz)});
        }
      }
    }
  };
  backup2D(specbleach2D_L, 0);
  backup2D(specbleach2D_R, 1);

  // Free existing handles
  if (specbleach1D_L)
    specbleach_free(specbleach1D_L);
  if (specbleach1D_R)
    specbleach_free(specbleach1D_R);
  if (specbleach2D_L)
    specbleach_2d_free(specbleach2D_L);
  if (specbleach2D_R)
    specbleach_2d_free(specbleach2D_R);

  // Initialize per-channel instances
  specbleach1D_L =
      specbleach_initialize(static_cast<uint32_t>(sampleRate), 50.0f);
  specbleach1D_R =
      specbleach_initialize(static_cast<uint32_t>(sampleRate), 50.0f);
  specbleach2D_L =
      specbleach_2d_initialize(static_cast<uint32_t>(sampleRate), 50.0f);
  specbleach2D_R =
      specbleach_2d_initialize(static_cast<uint32_t>(sampleRate), 50.0f);

  currentSampleRate = sampleRate;

  // Restore saved profiles into new handles
  for (const auto& item : profiles1D) {
    auto h = (item.channel == 1) ? specbleach1D_R : specbleach1D_L;
    if (h)
      loadProfile1DSafe(h, item.data.data(), item.size, item.blockCount,
                        item.mode);
  }

  for (const auto& item : profiles2D) {
    auto h = (item.channel == 1) ? specbleach2D_R : specbleach2D_L;
    if (h)
      loadProfile2DSafe(h, item.data.data(), item.size, item.blockCount,
                        item.mode);
  }

  // Restore state pending profiles if loaded before prepareToPlay
  bool hasCh1 = false;
  for (const auto& item : pendingProfiles) {
    if (item.channel == 1)
      hasCh1 = true;
    auto h1D = (item.channel == 1) ? specbleach1D_R : specbleach1D_L;
    auto h2D = (item.channel == 1) ? specbleach2D_R : specbleach2D_L;
    if (h1D)
      loadProfile1DSafe(h1D, item.data.data(), item.size, item.blockCount,
                        item.mode);
    if (h2D)
      loadProfile2DSafe(h2D, item.data.data(), item.size, item.blockCount,
                        item.mode);
  }

  // If no channel 1 profile, fallback channel 0 to channel 1
  if (!hasCh1) {
    for (const auto& item : pendingProfiles) {
      if (item.channel == 0) {
        if (specbleach1D_R)
          loadProfile1DSafe(specbleach1D_R, item.data.data(), item.size,
                            item.blockCount, item.mode);
        if (specbleach2D_R)
          loadProfile2DSafe(specbleach2D_R, item.data.data(), item.size,
                            item.blockCount, item.mode);
      }
    }
  }
  pendingProfiles.clear();

  // Synchronize loaded profiles between engines
  syncNoiseProfiles(0);
  syncNoiseProfiles(1);
}

void NoiseRepellentAudioProcessor::syncNoiseProfiles(int sourceAlgoMode) {
  auto syncHandle = [](auto getProfileFn, auto getSizeFn, auto getBlockCountFn,
                       auto isAvailableFn, auto loadFn, auto targetHandle,
                       auto sourceHandle) {
    if (!sourceHandle || !targetHandle)
      return;

    uint32_t sz = getSizeFn(sourceHandle);
    if (sz == 0)
      return;

    // Find a fallback profile (first available mode in source handle)
    float* fallbackP = nullptr;
    uint32_t fallbackBc = 0;

    for (int m = 1; m <= 4; ++m) {
      if (isAvailableFn(sourceHandle, m)) {
        fallbackP = getProfileFn(sourceHandle, m);
        fallbackBc = getBlockCountFn(sourceHandle, m);
        break;
      }
    }

    if (!fallbackP)
      return; // No profile learned in source handle yet

    // Load available mode profiles or use fallback so no target mode is left
    // empty
    for (int mode = 1; mode <= 4; ++mode) {
      if (isAvailableFn(sourceHandle, mode)) {
        float* p = getProfileFn(sourceHandle, mode);
        uint32_t bc = getBlockCountFn(sourceHandle, mode);
        if (p) {
          loadFn(targetHandle, p, sz, bc, mode);
        }
      } else {
        loadFn(targetHandle, fallbackP, sz, fallbackBc, mode);
      }
    }
  };

  if (sourceAlgoMode == 0) {
    // Sync 1D -> 2D
    syncHandle(specbleach_get_noise_profile_for_mode,
               specbleach_get_noise_profile_size,
               specbleach_get_noise_profile_block_count_for_mode,
               specbleach_noise_profile_available_for_mode,
               specbleach_2d_load_noise_profile_for_mode, specbleach2D_L,
               specbleach1D_L);

    syncHandle(specbleach_get_noise_profile_for_mode,
               specbleach_get_noise_profile_size,
               specbleach_get_noise_profile_block_count_for_mode,
               specbleach_noise_profile_available_for_mode,
               specbleach_2d_load_noise_profile_for_mode, specbleach2D_R,
               specbleach1D_R);
  } else if (sourceAlgoMode == 1) {
    // Sync 2D -> 1D
    syncHandle(specbleach_2d_get_noise_profile_for_mode,
               specbleach_2d_get_noise_profile_size,
               specbleach_2d_get_noise_profile_block_count_for_mode,
               specbleach_2d_noise_profile_available_for_mode,
               specbleach_load_noise_profile_for_mode, specbleach1D_L,
               specbleach2D_L);

    syncHandle(specbleach_2d_get_noise_profile_for_mode,
               specbleach_2d_get_noise_profile_size,
               specbleach_2d_get_noise_profile_block_count_for_mode,
               specbleach_2d_noise_profile_available_for_mode,
               specbleach_load_noise_profile_for_mode, specbleach1D_R,
               specbleach2D_R);
  }
}

void NoiseRepellentAudioProcessor::prepareToPlay(double sampleRate,
                                                 int samplesPerBlock) {
  ensureEnginesInitialized(sampleRate);

  // Set initial latency based on current algorithm mode
  currentAlgoMode = static_cast<int>(
      parameters.getRawParameterValue("algorithm_mode")->load());

  uint32_t latency = 0;
  if (currentAlgoMode == 0 && specbleach1D_L) {
    latency = specbleach_get_latency(specbleach1D_L);
  } else if (currentAlgoMode == 1 && specbleach2D_L) {
    latency = specbleach_2d_get_latency(specbleach2D_L);
  }
  setLatencySamples(static_cast<int>(latency));

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
  dryWetMixer.setWetLatency(static_cast<float>(latency));

  currentLatency = latency;
  visualizerDelayBuffer.assign(std::max<size_t>(32768, latency + 8192), 0.0f);
  visualizerDelayWritePos = 0;

  // Pre-allocate persistent buffers to prevent audio-thread allocations
  dryInputL.resize(static_cast<size_t>(bufferCapacity), 0.0f);
  crossfadeBuffer.setSize(static_cast<int>(preparedNumChannels), bufferCapacity,
                          false, false, true);
  crossfadeDelayBuffer.setSize(static_cast<int>(preparedNumChannels), 16384,
                               false, false, true);
  crossfadeDelayBuffer.clear();
  crossfadeDelayWritePos = 0;
  crossfadeLatencyDiff = 0;
  delaySlewProgress = 1.0f;
  delaySlewStep = 0.0f;
  jassert(dryInputL.size() >= static_cast<size_t>(samplesPerBlock));
  jassert(crossfadeBuffer.getNumSamples() >= samplesPerBlock);
  crossfadeProgress = 1.0f;
  sourceAlgoMode = currentAlgoMode;
  targetAlgoMode = currentAlgoMode;

  // Reset FFT accumulation
  fftAccumInput.fill(0.0f);
  fftAccumOutput.fill(0.0f);
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

  const int safeNumChannels =
      std::min(numChannels, static_cast<int>(preparedNumChannels));
  const int safeNumSamples =
      std::min(numSamples, static_cast<int>(preparedBlockSize));

  const bool isBypassed =
      parameters.getRawParameterValue("bypass")->load() > 0.5f;
  const int algoMode = static_cast<int>(
      parameters.getRawParameterValue("algorithm_mode")->load());
  const int hpssQuality =
      static_cast<int>(parameters.getRawParameterValue("hpss_quality")->load());
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

  // Suppression: pass raw 0-100 value — libspecbleach divides by 100 internally
  const float suppression =
      parameters.getRawParameterValue("suppression_strength")->load();

  const bool residualListen =
      parameters.getRawParameterValue("residual_listen")->load() > 0.5f;
  const float profileOffset =
      parameters.getRawParameterValue("noise_profile_offset")->load();
  const bool curveEnabled =
      parameters.getRawParameterValue("reduction_curve_enabled")->load() > 0.5f;

  if (curveEnabled) {
    if (curveNodesDirty.exchange(false)) {
      juce::SpinLock::ScopedTryLockType tryLock(curveLock);
      if (tryLock.isLocked()) {
        uint32_t profileSize = 0;
        if (algoMode == 0 && specbleach1D_L)
          profileSize = specbleach_get_noise_profile_size(specbleach1D_L);
        else if (algoMode == 1 && specbleach2D_L)
          profileSize = specbleach_2d_get_noise_profile_size(specbleach2D_L);
        if (profileSize > 0)
          interpolateCurve(profileSize);
      } else {
        curveNodesDirty = true;
      }
    }
  }

  // Sync profiles if manual noise learning was just turned off
  if (wasLearning && !learnNoise) {
    juce::SpinLock::ScopedTryLockType tryLock(profileLock);
    if (tryLock.isLocked()) {
      syncNoiseProfiles(currentAlgoMode);
    }
  }
  wasLearning = learnNoise;

  // Prepare parameters for DSP engines
  SpectralBleachDenoiserParameters p;
  p.learn_noise = learnNoise ? 1 : 0;
  p.residual_listen = residualListen;
  p.reduction_amount = masterRed;
  p.smoothing_factor = smoothing;
  p.whitening_factor = whitening;
  p.adaptive_noise = adaptiveNoise ? 1 : 0;
  p.noise_estimation_method = adaptiveMethod;
  p.masking_depth = masking;
  p.suppression_strength = suppression;
  p.aggressiveness = aggressiveness;
  p.tonal_reduction = tonalRed;
  p.hpss_quality_mode = hpssQuality;
  p.noise_profile_offset_db = profileOffset;
  p.reduction_curve_enabled = curveEnabled;
  p.reduction_curve_bias =
      curveEnabled ? interpolatedCurveBias.data() : nullptr;

  SpectralBleach2DDenoiserParameters p2;
  p2.learn_noise = learnNoise ? 1 : 0;
  p2.residual_listen = residualListen;
  p2.reduction_amount = masterRed;
  p2.smoothing_factor = smoothing;
  p2.whitening_factor = whitening;
  p2.adaptive_noise = adaptiveNoise ? 1 : 0;
  p2.noise_estimation_method = adaptiveMethod;
  p2.nlm_masking_protection = masking;
  p2.suppression_strength = suppression;
  p2.aggressiveness = aggressiveness;
  p2.tonal_reduction = tonalRed;
  p2.hpss_quality_mode = hpssQuality;
  p2.noise_profile_offset_db = profileOffset;
  p2.reduction_curve_enabled = curveEnabled;
  p2.reduction_curve_bias =
      curveEnabled ? interpolatedCurveBias.data() : nullptr;

  // Update latency dynamically, sync profiles, and start smooth crossfade when
  // algorithm mode changes
  if (algoMode != currentAlgoMode) {
    juce::SpinLock::ScopedTryLockType tryLock(profileLock);
    if (tryLock.isLocked()) {
      syncNoiseProfiles(currentAlgoMode);
    }
    sourceAlgoMode = currentAlgoMode;
    targetAlgoMode = algoMode;
    currentAlgoMode = algoMode;

    if (specbleach1D_L)
      specbleach_load_parameters(specbleach1D_L, p);
    if (specbleach1D_R)
      specbleach_load_parameters(specbleach1D_R, p);
    if (specbleach2D_L)
      specbleach_2d_load_parameters(specbleach2D_L, p2);
    if (specbleach2D_R)
      specbleach_2d_load_parameters(specbleach2D_R, p2);

    uint32_t lat1D =
        specbleach1D_L ? specbleach_get_latency(specbleach1D_L) : 0;
    uint32_t lat2D =
        specbleach2D_L ? specbleach_2d_get_latency(specbleach2D_L) : 0;
    crossfadeLatencyDiff = (lat2D > lat1D) ? (lat2D - lat1D) : 0;

    // Initiate 40 ms smooth crossfade transition to prevent clicks/pops
    int crossfadeSamples = static_cast<int>(currentSampleRate * 0.040);
    if (crossfadeSamples < 1)
      crossfadeSamples = 1;
    crossfadeProgress = 0.0f;
    crossfadeStep = 1.0f / static_cast<float>(crossfadeSamples);

    if (targetAlgoMode == 0) {
      delaySlewProgress = 0.0f;
      delaySlewStep = 1.0f / static_cast<float>(crossfadeSamples);
    } else {
      delaySlewProgress = 1.0f;
      delaySlewStep = 0.0f;
    }
  }

  if (hpssQuality != currentHpssQuality) {
    currentHpssQuality = hpssQuality;
    if (specbleach1D_L)
      specbleach_load_parameters(specbleach1D_L, p);
    if (specbleach1D_R)
      specbleach_load_parameters(specbleach1D_R, p);
    if (specbleach2D_L)
      specbleach_2d_load_parameters(specbleach2D_L, p2);
    if (specbleach2D_R)
      specbleach_2d_load_parameters(specbleach2D_R, p2);
  }

  // Save dry input copy for FFT visualization before processing
  const size_t copySamples =
      std::min(static_cast<size_t>(numSamples), dryInputL.size());
  if (numChannels >= 1 && copySamples > 0)
    std::copy_n(buffer.getReadPointer(0), copySamples, dryInputL.begin());

  // Push dry samples into JUCE DryWetMixer before in-place DSP processing
  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);

  if (crossfadeProgress < 1.0f) {
    // Smooth crossfade mode: run both engines and blend sample-by-sample with
    // delay alignment
    jassert(numSamples <= crossfadeBuffer.getNumSamples());
    const int maxChannels =
        std::min(numChannels, crossfadeBuffer.getNumChannels());

    for (int ch = 0; ch < maxChannels; ++ch)
      std::copy_n(buffer.getReadPointer(ch), numSamples,
                  crossfadeBuffer.getWritePointer(ch));

    // Process 1D on buffer
    if (specbleach1D_L) {
      specbleach_load_parameters(specbleach1D_L, p);
      float* ch0 = buffer.getWritePointer(0);
      specbleach_process(specbleach1D_L, static_cast<uint32_t>(numSamples), ch0,
                         ch0);
    }
    if (numChannels >= 2 && specbleach1D_R) {
      specbleach_load_parameters(specbleach1D_R, p);
      float* ch1 = buffer.getWritePointer(1);
      specbleach_process(specbleach1D_R, static_cast<uint32_t>(numSamples), ch1,
                         ch1);
    }

    // Process 2D on crossfadeBuffer
    if (maxChannels >= 1 && specbleach2D_L) {
      specbleach_2d_load_parameters(specbleach2D_L, p2);
      float* ch0 = crossfadeBuffer.getWritePointer(0);
      specbleach_2d_process(specbleach2D_L, static_cast<uint32_t>(numSamples),
                            ch0, ch0);
    }
    if (numChannels >= 2 && specbleach2D_R && maxChannels >= 2) {
      specbleach_2d_load_parameters(specbleach2D_R, p2);
      float* ch1 = crossfadeBuffer.getWritePointer(1);
      specbleach_2d_process(specbleach2D_R, static_cast<uint32_t>(numSamples),
                            ch1, ch1);
    }

    // Equal-power crossfade blend with delay alignment:
    // Delay 1D output by crossfadeLatencyDiff so it precisely matches 2D
    // latency during crossfade
    const size_t delayBufSize =
        static_cast<size_t>(crossfadeDelayBuffer.getNumSamples());
    size_t writePos = crossfadeDelayWritePos;

    for (int ch = 0; ch < maxChannels; ++ch) {
      float* out1D = buffer.getWritePointer(ch);
      const float* out2D = crossfadeBuffer.getReadPointer(ch);
      float* delayLine = crossfadeDelayBuffer.getWritePointer(ch);

      float currentP = crossfadeProgress;
      size_t wp = writePos;

      for (int s = 0; s < numSamples; ++s) {
        float in1D = out1D[s];
        delayLine[wp] = in1D;
        size_t delay = crossfadeLatencyDiff % delayBufSize;
        size_t readPos = (wp + delayBufSize - delay) % delayBufSize;
        float aligned1D = delayLine[readPos];
        wp = (wp + 1) % delayBufSize;

        float pVal = std::min(1.0f, currentP);
        float rad = pVal * juce::MathConstants<float>::halfPi;
        // If transitioning 1D -> 2D (target == 1): fade from aligned 1D to 2D
        // If transitioning 2D -> 1D (target == 0): fade from 2D to aligned 1D
        float w1D = (targetAlgoMode == 1) ? std::cos(rad) : std::sin(rad);
        float w2D = (targetAlgoMode == 1) ? std::sin(rad) : std::cos(rad);

        out1D[s] = w1D * aligned1D + w2D * out2D[s];
        currentP += crossfadeStep;
      }
    }
    crossfadeDelayWritePos = (writePos + numSamples) % delayBufSize;
    crossfadeProgress += numSamples * crossfadeStep;
    if (crossfadeProgress > 1.0f)
      crossfadeProgress = 1.0f;
  } else if (algoMode == 0) {
    // 1D Mode
    if (specbleach1D_L) {
      specbleach_load_parameters(specbleach1D_L, p);
      float* ch0 = buffer.getWritePointer(0);
      specbleach_process(specbleach1D_L, static_cast<uint32_t>(numSamples), ch0,
                         ch0);
    }
    if (numChannels >= 2 && specbleach1D_R) {
      specbleach_load_parameters(specbleach1D_R, p);
      float* ch1 = buffer.getWritePointer(1);
      specbleach_process(specbleach1D_R, static_cast<uint32_t>(numSamples), ch1,
                         ch1);
    }

    if (delaySlewProgress < 1.0f) {
      const size_t delayBufSize =
          static_cast<size_t>(crossfadeDelayBuffer.getNumSamples());
      const int maxChannels =
          std::min(numChannels, crossfadeDelayBuffer.getNumChannels());
      size_t writePos = crossfadeDelayWritePos;

      for (int ch = 0; ch < maxChannels; ++ch) {
        float* out1D = buffer.getWritePointer(ch);
        float* delayLine = crossfadeDelayBuffer.getWritePointer(ch);
        float currentP = delaySlewProgress;
        size_t wp = writePos;

        for (int s = 0; s < numSamples; ++s) {
          float undelayedSample = out1D[s];
          delayLine[wp] = undelayedSample;
          size_t delay = crossfadeLatencyDiff % delayBufSize;
          size_t readPos = (wp + delayBufSize - delay) % delayBufSize;
          float alignedSample = delayLine[readPos];
          wp = (wp + 1) % delayBufSize;

          float pVal = std::min(1.0f, currentP);
          float rad = pVal * juce::MathConstants<float>::halfPi;
          // Smoothly fade from aligned (delayed) to undelayed 1D
          float wDelayed = std::cos(rad);
          float wUndelayed = std::sin(rad);

          out1D[s] = wDelayed * alignedSample + wUndelayed * undelayedSample;
          currentP += delaySlewStep;
        }
      }
      crossfadeDelayWritePos = (writePos + numSamples) % delayBufSize;
      delaySlewProgress += numSamples * delaySlewStep;
      if (delaySlewProgress > 1.0f)
        delaySlewProgress = 1.0f;
    }
  } else if (algoMode == 1) {
    // Steady-state 2D mode
    if (specbleach2D_L) {
      specbleach_2d_load_parameters(specbleach2D_L, p2);
      float* ch0 = buffer.getWritePointer(0);
      specbleach_2d_process(specbleach2D_L, static_cast<uint32_t>(numSamples),
                            ch0, ch0);
    }
    if (numChannels >= 2 && specbleach2D_R) {
      specbleach_2d_load_parameters(specbleach2D_R, p2);
      float* ch1 = buffer.getWritePointer(1);
      specbleach_2d_process(specbleach2D_R, static_cast<uint32_t>(numSamples),
                            ch1, ch1);
    }
  }

  // Update dynamic latency if changed by algorithm mode or HPSS quality mode
  uint32_t activeLatency = 0;
  if (algoMode == 0 && specbleach1D_L) {
    activeLatency = specbleach_get_latency(specbleach1D_L);
  } else if (algoMode == 1 && specbleach2D_L) {
    activeLatency = specbleach_2d_get_latency(specbleach2D_L);
  }

  if (activeLatency != currentLatency) {
    currentLatency = activeLatency;
    setLatencySamples(static_cast<int>(activeLatency));
    dryWetMixer.setWetLatency(static_cast<float>(activeLatency));
    if (visualizerDelayBuffer.size() < activeLatency + 8192) {
      visualizerDelayBuffer.resize(activeLatency + 8192, 0.0f);
    }
  }

  // Apply Soft Crossfade Bypass using JUCE DryWetMixer (with dry latency
  // compensation)
  dryWetMixer.setWetMixProportion(isBypassed ? 0.0f : 1.0f);
  dryWetMixer.mixWetSamples(audioBlock);

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

    if (fftAccumCount < kFftSize) {
      fftAccumInput[fftAccumCount] = alignedInputSample;
      fftAccumOutput[fftAccumCount] = outputSrc ? outputSrc[s] : 0.0f;
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
        bool profileHasAnyMode = false;

        juce::SpinLock::ScopedTryLockType profileTryLock(profileLock);
        if (profileTryLock.isLocked()) {
          if (algoMode == 0 && specbleach1D_L) {
            for (int m = 1; m <= 4; ++m) {
              if (specbleach_noise_profile_available_for_mode(specbleach1D_L,
                                                              m)) {
                profileHasAnyMode = true;
                break;
              }
            }
            if (isAdaptive || profileHasAnyMode) {
              actualNoiseProfile =
                  specbleach_get_active_noise_profile(specbleach1D_L);
              profileSize = specbleach_get_noise_profile_size(specbleach1D_L);
              if (!actualNoiseProfile && profileHasAnyMode) {
                for (int m = 1; m <= 4; ++m) {
                  if (specbleach_noise_profile_available_for_mode(
                          specbleach1D_L, m)) {
                    actualNoiseProfile = specbleach_get_noise_profile_for_mode(
                        specbleach1D_L, m);
                    break;
                  }
                }
              }
              profileAvailable =
                  (actualNoiseProfile != nullptr && profileSize > 0);
            }
          } else if (algoMode == 1 && specbleach2D_L) {
            for (int m = 1; m <= 4; ++m) {
              if (specbleach_2d_noise_profile_available_for_mode(specbleach2D_L,
                                                                 m)) {
                profileHasAnyMode = true;
                break;
              }
            }
            if (isAdaptive || profileHasAnyMode) {
              actualNoiseProfile =
                  specbleach_2d_get_active_noise_profile(specbleach2D_L);
              profileSize =
                  specbleach_2d_get_noise_profile_size(specbleach2D_L);
              if (!actualNoiseProfile && profileHasAnyMode) {
                for (int m = 1; m <= 4; ++m) {
                  if (specbleach_2d_noise_profile_available_for_mode(
                          specbleach2D_L, m)) {
                    actualNoiseProfile =
                        specbleach_2d_get_noise_profile_for_mode(specbleach2D_L,
                                                                 m);
                    break;
                  }
                }
              }
              profileAvailable =
                  (actualNoiseProfile != nullptr && profileSize > 0);
            }
          }
        }

        frame.hasNoiseProfile = profileAvailable;
        frame.isLinked =
            (parameters.getRawParameterValue("link_reduction")->load() > 0.5f);
        frame.reductionCurveEnabled =
            (parameters.getRawParameterValue("reduction_curve_enabled")
                 ->load() > 0.5f);
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
          std::array<float, 16> peakBuf{};
          uint32_t numPeaks = 0;
          if (algoMode == 0 && specbleach1D_L) {
            numPeaks = specbleach_get_tonal_peaks_for_profile(
                specbleach1D_L, actualNoiseProfile,
                static_cast<uint32_t>(realProfileBins), peakBuf.data(),
                static_cast<uint32_t>(peakBuf.size()));
          } else if (algoMode == 1 && specbleach2D_L) {
            numPeaks = specbleach_2d_get_tonal_peaks_for_profile(
                specbleach2D_L, actualNoiseProfile,
                static_cast<uint32_t>(realProfileBins), peakBuf.data(),
                static_cast<uint32_t>(peakBuf.size()));
          }

          for (uint32_t i = 0; i < numPeaks; ++i) {
            frame.tonalPeaksHz.push_back(peakBuf[i]);
          }
        } else {
          frame.noiseFloorDB.fill(-120.0f);
        }

        spectralFifo.finishedWrite(1);
      }

      // Shift: keep last quarter for overlap (hop = 75%)
      constexpr size_t kHop = kFftSize / 4;
      std::memmove(fftAccumInput.data(), fftAccumInput.data() + kHop,
                   (kFftSize - kHop) * sizeof(float));
      std::memmove(fftAccumOutput.data(), fftAccumOutput.data() + kHop,
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
  if (specbleach1D_L)
    specbleach_reset_noise_profile(specbleach1D_L);
  if (specbleach1D_R)
    specbleach_reset_noise_profile(specbleach1D_R);
  if (specbleach2D_L)
    specbleach_2d_reset_noise_profile(specbleach2D_L);
  if (specbleach2D_R)
    specbleach_2d_reset_noise_profile(specbleach2D_R);
  pendingProfiles.clear();
}

bool NoiseRepellentAudioProcessor::hasNoiseProfile() const {
  for (int mode = 1; mode <= 4; ++mode) {
    if (specbleach1D_L &&
        specbleach_noise_profile_available_for_mode(specbleach1D_L, mode))
      return true;
    if (specbleach2D_L &&
        specbleach_2d_noise_profile_available_for_mode(specbleach2D_L, mode))
      return true;
  }
  return false;
}

juce::AudioProcessorEditor* NoiseRepellentAudioProcessor::createEditor() {
  return new NoiseRepellentAudioProcessorEditor(*this);
}

void NoiseRepellentAudioProcessor::getStateInformation(
    juce::MemoryBlock& destData) {
  juce::SpinLock::ScopedLockType lock(profileLock);

  // Sync profiles across engines before state serialization
  syncNoiseProfiles(currentAlgoMode);

  auto state = parameters.copyState();

  juce::ValueTree profilesTree("LEARNED_PROFILES");

  for (int channel = 0; channel < 2; ++channel) {
    auto h1D = (channel == 0) ? specbleach1D_L : specbleach1D_R;
    auto h2D = (channel == 0) ? specbleach2D_L : specbleach2D_R;

    for (int mode = 1; mode <= 4; ++mode) {
      const float* profile = nullptr;
      uint32_t profileSize = 0;
      uint32_t blockCount = 0;

      if (h1D && specbleach_noise_profile_available_for_mode(h1D, mode)) {
        profile = specbleach_get_noise_profile_for_mode(h1D, mode);
        profileSize = specbleach_get_noise_profile_size(h1D);
        blockCount =
            specbleach_get_noise_profile_block_count_for_mode(h1D, mode);
      } else if (h2D &&
                 specbleach_2d_noise_profile_available_for_mode(h2D, mode)) {
        profile = specbleach_2d_get_noise_profile_for_mode(h2D, mode);
        profileSize = specbleach_2d_get_noise_profile_size(h2D);
        blockCount =
            specbleach_2d_get_noise_profile_block_count_for_mode(h2D, mode);
      }

      if (profile && profileSize > 0) {
        juce::MemoryBlock mb(profile, profileSize * sizeof(float));
        juce::ValueTree node("CHANNEL_PROFILE");
        node.setProperty("channel", channel, nullptr);
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

            auto h1D = (channel == 1) ? specbleach1D_R : specbleach1D_L;
            auto h2D = (channel == 1) ? specbleach2D_R : specbleach2D_L;

            if (h1D)
              loadProfile1DSafe(h1D, floatArray, size, blockCount, mode);
            if (h2D)
              loadProfile2DSafe(h2D, floatArray, size, blockCount, mode);
          }
        }
      }

      // If legacy profile had no Channel 1 profile, copy Channel 0 to Channel 1
      // handles
      if (!hasCh1) {
        for (const auto& pp : pendingProfiles) {
          if (pp.channel == 0) {
            if (specbleach1D_R)
              loadProfile1DSafe(specbleach1D_R, pp.data.data(), pp.size,
                                pp.blockCount, pp.mode);
            if (specbleach2D_R)
              loadProfile2DSafe(specbleach2D_R, pp.data.data(), pp.size,
                                pp.blockCount, pp.mode);
          }
        }
      }

      // Sync loaded profiles to ensure both engines are fully populated
      syncNoiseProfiles(0);
      syncNoiseProfiles(1);

      if (specbleach1D_L != nullptr) {
        pendingProfiles.clear();
      }
    }
  }
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
  return new NoiseRepellentAudioProcessor();
}
