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
#include "specbleach_version.h"
#include <algorithm>
#include <cstddef>
#include <memory>

namespace {

// Silence gate: amplitude floor ~-100 dB, i.e. the amplitude-domain twin of
// libspecbleach's ESTIMATOR_SILENCE_THRESHOLD (1e-10 power).
constexpr float kSilenceAmplitudeFloor = 1e-5f;

// ---------------------------------------------------------------------------
// File-local helpers for the parameter layout, state parsing and the
// processBlock gating. Free functions keep PluginProcessor.h minimal; the
// member state they need is passed by reference.
// ---------------------------------------------------------------------------

juce::AudioProcessorValueTreeState::ParameterLayout
createParameterLayout() {
  std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

  params.push_back(std::make_unique<juce::AudioParameterChoice>(
      "algorithm_mode", "Smoothing Method",
      juce::StringArray{"Temporal (Fast / Low CPU)",
                        "Patch-Based NLM (Texture & Detail)",
                        "NLM + DFTT Refinement (Tonal / Difficult)",
                        "Collaborative BM3D (Smooth / Steady)"},
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

// Fallback state parser for hosts that hand back plain XML text instead of the
// binary produced by copyXmlToBinary.
std::unique_ptr<juce::XmlElement> parsePlainXml(const std::byte* data,
                                                int sizeInBytes) {
  if (data == nullptr || sizeInBytes <= 0)
    return nullptr;
  const juce::String xmlString =
      juce::String::createStringFromData(data, sizeInBytes);
  if (xmlString.trim().startsWith("<"))
    return juce::XmlDocument::parse(xmlString);
  return nullptr;
}

// Transport state for the GUI; unknown when the host provides no position.
void updateTransportFromHost(juce::AudioPlayHead* playHead,
                             std::atomic<bool>& playing,
                             std::atomic<bool>& stateKnown) {
  if (playHead == nullptr) {
    stateKnown.store(false);
    return;
  }
  if (const auto position = playHead->getPosition()) {
    playing.store(position->getIsPlaying());
    stateKnown.store(true);
  } else {
    stateKnown.store(false);
  }
}

// Silence detection: when input is ~-100dB and not learning, skip heavy
// STFT+denoise (was culprit for idle > denoising and no-playback CPU).
bool detectSilence(const juce::AudioBuffer<float>& buffer, bool learnNoise,
                   bool adaptiveNoise) {
  if (learnNoise)
    return false;
  const float thresh = kSilenceAmplitudeFloor;
  bool isSilent = true;
  const int numChannels = buffer.getNumChannels();
  const int numSamples = buffer.getNumSamples();
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
  return isSilent && !adaptiveNoise;
}

// Advances the per-block silence streak and reports whether the pipeline is
// fully flushed (reported latency plus gain-release tails) so the engine may
// sleep; sleeping earlier would chop decaying tails.
bool engineCanSleep(bool isSilent, bool hasProfile, bool learnNoise,
                    bool adaptiveNoise, uint32_t& silentBlocksStreak,
                    uint32_t latency, double sampleRate, int numSamples) {
  if (isSilent)
    ++silentBlocksStreak;
  else
    silentBlocksStreak = 0;
  const double tailMarginSamples =
      (sampleRate > 0.0 ? sampleRate * 0.5 : 24000.0);
  const double flushSamples = static_cast<double>(latency) + tailMarginSamples;
  const bool pipelineFlushed =
      static_cast<double>(silentBlocksStreak) *
          static_cast<double>(std::max(numSamples, 1)) >
      flushSamples;
  return isSilent && pipelineFlushed && !hasProfile && !adaptiveNoise &&
         !learnNoise;
}

} // namespace

NoiseRepellentAudioProcessor::NoiseRepellentAudioProcessor()
    : AudioProcessor(
          BusesProperties()
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      parameters(*this, nullptr, "PARAMETERS", createParameterLayout()) {
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
  dryWetMixer.reset();
}

juce::AudioProcessorParameter*
NoiseRepellentAudioProcessor::getBypassParameter() const {
  if (wrapperType == wrapperType_LV2) {
    // JUCE's LV2 wrapper mirrors lv2:enabled into this parameter on every run.
    // Keep the plugin's own bypass control independent; the wrapper calls
    // processBlockBypassed() for host bypass when no parameter is returned.
    return nullptr;
  }
  return bypassParameter;
}

float NoiseRepellentAudioProcessor::getFrameSizeMs() const {
  if (const auto* choice = dynamic_cast<juce::AudioParameterChoice*>(
          parameters.getParameter("frame_size"))) {
    const int index = choice->getIndex();
    if (index >= 0 && index < 5)
      return kFrameSizeOptionsMs[static_cast<size_t>(index)];
  }
  return kFrameSizeOptionsMs[kDefaultFrameSizeIndex];
}

bool NoiseRepellentAudioProcessor::isLowLatency() const {
  if (const auto* p = parameters.getRawParameterValue("low_latency"))
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
    // ponytail: coalesce — LV2 delivers structural-param automation on the
    // audio thread per sub-block; one flag beats a callAsync storm that
    // resets the engine under an in-flight run().
    // A change driven by a state restore (GH #221) must keep the profiles
    // the restore is about to install when this deferred rebuild fires.
    if (status.replacingStateForRestore.load())
      status.preserveProfilesOnFrameRebuild.store(true);
    status.frameSizeRebuildPending.store(true);
  }
}

void NoiseRepellentAudioProcessor::rebuildForFrameSizeChange(
    bool onAudioThread) {
  if (status.currentSampleRate <= 0.0 || engine.group() == nullptr)
    return; // not prepared yet — next prepareToPlay picks up the new value
  // Auto-stop an in-progress Learn: a half-rolled mean must not migrate
  // across resolutions.
  if (auto* learn = parameters.getParameter("learn_noise"))
    learn->setValueNotifyingHost(0.0f);
  // Low-latency transitions install control defaults: on entry smoothing
  // 30, full aggressiveness, no masking veto; on exit the parameter
  // defaults return. Applied once per transition (currentLowLatency still
  // holds the old value here); the user stays free to adjust afterwards.
  auto pushValue = [this](const char* id, float denormalized) {
    if (auto* p = parameters.getParameter(id))
      p->setValueNotifyingHost(p->convertTo0to1(denormalized));
  };
  auto pushDefault = [this](const char* id) {
    if (auto* p = parameters.getParameter(id))
      p->setValueNotifyingHost(p->getDefaultValue());
  };
  if (isLowLatency() && !status.currentLowLatency) {
    pushValue("smoothing_factor", 30.0f);
    pushValue("aggressiveness", 1.0f);
    pushValue("masking_depth", 0.0f);
  } else if (!isLowLatency() && status.currentLowLatency) {
    pushDefault("smoothing_factor");
    pushDefault("aggressiveness");
    pushDefault("masking_depth");
  }
  // On the audio thread the rebuild is serialized with the engine run in
  // processBlock, so no suspend dance is needed (that only guards
  // message-thread entry).
  if (!onAudioThread)
    suspendProcessing(true);
  const uint32_t wantedChannels = static_cast<uint32_t>(
      std::max({getTotalNumInputChannels(), getTotalNumOutputChannels(), 1}));
  engine.ensureInitialized(status.currentSampleRate, wantedChannels,
                           isLowLatency(), getFrameSizeMs());
  updateLatencyReporting();
  // Fresh STFT history: clear visualization accumulators so the first frames
  // after the switch don't mix pre/post-switch spectra. The rebuilt
  // pipeline is empty, so the silence streak restarts clean.
  viz.resetAccumulators();
  status.silentBlocksStreak = 0;
  if (!onAudioThread)
    suspendProcessing(false);
  if (auto* editor = getActiveEditor())
    editor->repaint();
}

void NoiseRepellentAudioProcessor::updateLatencyReporting() {
  // Latency depends on the STFT frame size (frame + NLM look-ahead) but is
  // constant across smoothing modes (temporal is padded to the NLM
  // look-ahead), so hosts only need a re-report after a frame-size switch.
  // Low-latency mode is the second exception: causal 1D-only, zero
  // look-ahead (latency = frame alone), rebuilt via the same suspend path.
  const uint32_t reportedLatency =
      engine.group() != nullptr ? specbleach_stereo_get_latency(engine.group())
                                : 0;
  status.lastReportedLatency = reportedLatency;
  setLatencySamples(static_cast<int>(reportedLatency));
  dryWetMixer.setWetLatency(static_cast<float>(reportedLatency));

  status.currentLatency = reportedLatency;
  viz.resizeDelayLine(reportedLatency);
}

void NoiseRepellentAudioProcessor::prepareToPlay(double sampleRate,
                                                 int samplesPerBlock) {
  const uint32_t wantedChannels = static_cast<uint32_t>(
      std::max({getTotalNumInputChannels(), getTotalNumOutputChannels(), 1}));
  engine.ensureInitialized(sampleRate, wantedChannels, isLowLatency(),
                           getFrameSizeMs());

  const int bufferCapacity = std::max(samplesPerBlock, 16384);
  status.preparedBlockSize = static_cast<uint32_t>(bufferCapacity);
  status.preparedNumChannels = wantedChannels;

  juce::dsp::ProcessSpec spec;
  spec.sampleRate = sampleRate;
  spec.maximumBlockSize = status.preparedBlockSize;
  spec.numChannels = status.preparedNumChannels;

  dryWetMixer.prepare(spec);
  dryWetMixer.setMixingRule(juce::dsp::DryWetMixingRule::linear);
  updateLatencyReporting();

  // Pre-allocate persistent buffers to prevent audio-thread allocations
  dryInputL.resize(static_cast<size_t>(bufferCapacity), 0.0f);
  jassert(dryInputL.size() >= static_cast<size_t>(samplesPerBlock));

  viz.prepareFrames();

  // NB: pendingProfiles is owned by EngineRuntime::ensureInitialized —
  // retained (narrow-group) entries must survive prepareToPlay for a later
  // rebuild.

  viz.resetAccumulators();
  // Unknown pipeline state after (re)configuration: restart the silence
  // streak so a fresh silence re-proves the flush before any sleep.
  status.silentBlocksStreak = 0;
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
  updateTransportFromHost(getPlayHead(), status.hostTransportPlaying,
                          status.hostTransportStateKnown);

  const int numSamples = buffer.getNumSamples();
  const int numChannels = buffer.getNumChannels();

  if (numSamples == 0 || numChannels == 0)
    return;

  // Track genuinely in-progress offline renders (see getNonRealtimeBlockCount).
  if (isNonRealtime())
    status.nonRealtimeBlockCount.fetch_add(1);

  // Single coalesced structural rebuild, serialized with the engine run
  // below — the race-free funnel for off-message-thread
  // frame_size/low_latency changes (LV2 automation). Never concurrent, never
  // a storm.
  if (status.frameSizeRebuildPending.exchange(false))
    rebuildForFrameSizeChange(/*onAudioThread=*/true);

  const bool isBypassed =
      parameters.getRawParameterValue("bypass")->load() > 0.5f;
  const EngineRuntime::EngineParams ep = engine.buildParams(parameters);

  // Silence detection: when input is ~-100dB and not learning, skip heavy
  // STFT+denoise (was culprit for idle > denoising and no-playback CPU).
  const bool isSilent = detectSilence(buffer, ep.learnNoise, ep.adaptiveNoise);

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
  if (engineCanSleep(isSilent, hasProfile, ep.learnNoise, ep.adaptiveNoise,
                     status.silentBlocksStreak, status.currentLatency,
                     status.currentSampleRate, numSamples) &&
      !isBypassed) {
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
    engine.run(buffer, ep);
  }

  // Query transient protection status and intensity (aggregated max across
  // channels)
  const float reportedTransientIntensity = engine.transientIntensity();
  status.transientActivity.store(reportedTransientIntensity);
  status.transientProtectionActive.store(ep.transientProtectionEnable);

  // Apply Soft Crossfade Bypass using JUCE DryWetMixer (with dry latency
  // compensation)
  dryWetMixer.setWetMixProportion(isBypassed ? 0.0f : 1.0f);
  dryWetMixer.mixWetSamples(audioBlock);

  // Skip FFT analysis and FIFO writes during offline rendering or if the GUI
  // is closed
  if (isNonRealtime() || getActiveEditor() == nullptr) {
    viz.resetAccumCount();
    return;
  }

  VisualizationFeed::FeedInput vizInput{
      .group = engine.group(),
      .profileLock = profiles.lock(),
      .apvts = parameters,
      .latency = status.currentLatency,
      .isLearning = ep.learnNoise,
      .adaptiveNoise = ep.adaptiveNoise,
      .transientIntensity = reportedTransientIntensity,
      .transientProtectionActive = ep.transientProtectionEnable};

  // When input is silent and not learning, skip heavy 4096 FFT but push a
  // cheap silent frame so input/output fall to -120 dB while the learned
  // noise profile stays frozen (only reset clears it).
  if (isSilent) {
    viz.resetAccumCount();
    viz.pushSilentFrame(vizInput, profiles);
    return;
  }

  // ── FFT-based spectral frame for GUI visualization ──
  // Accumulate samples until we have a full FFT window (kFftSize samples)
  viz.feed(vizInput, dryInputL.data(),
           (numChannels >= 1) ? buffer.getReadPointer(0) : nullptr,
           copySamples);
}

void NoiseRepellentAudioProcessor::processBlockBypassed(
    juce::AudioBuffer<float>& buffer, juce::MidiBuffer& /*midiMessages*/) {
  if (const int numSamples = buffer.getNumSamples(),
      numChannels = buffer.getNumChannels();
      numSamples == 0 || numChannels == 0)
    return;

  if (isNonRealtime())
    status.nonRealtimeBlockCount.fetch_add(1);

  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);
  // Run the engine even though the wet path is muted: same alignment
  // requirement as the internal bypass (see processBlock) — the buffer
  // must hold pipeline-delayed content while the mixer ramps down, and
  // the delay line stays warm for the resume boundary.
  engine.run(buffer, engine.buildParams(parameters));
  dryWetMixer.setWetMixProportion(0.0f);
  dryWetMixer.mixWetSamples(audioBlock);
}

void NoiseRepellentAudioProcessor::getStateInformation(
    juce::MemoryBlock& destData) {
  juce::SpinLock::ScopedLockType lock(profiles.lock());

  auto state = parameters.copyState();

  // Fill per-channel/per-mode gaps before serialization (inside
  // appendToState), then the REDUCTION_CURVE nodes.
  profiles.appendToState(engine.group(), state);
  curve.appendToState(state);

  std::unique_ptr<juce::XmlElement> xml(state.createXml());
  copyXmlToBinary(*xml, destData);
}

void NoiseRepellentAudioProcessor::setStateInformation(const void* data,
                                                       int sizeInBytes) {
  // The JUCE override signature mandates void*; hand the helpers a
  // meaningful byte-pointer type instead (cpp:S5008).
  const auto* stateBytes = static_cast<const std::byte*>(data);
  profiles.clear();

  std::unique_ptr<juce::XmlElement> xmlState(
      getXmlFromBinary(data, sizeInBytes));
  if (xmlState == nullptr)
    xmlState = parsePlainXml(stateBytes, sizeInBytes);

  if (xmlState != nullptr) {
    juce::ValueTree state = juce::ValueTree::fromXml(*xmlState);

    // Extract REDUCTION_CURVE if present
    std::vector<CurveNode> loadedNodes;
    if (ReductionCurveState::decodeNodesFrom(state, *xmlState, loadedNodes))
      curve.setNodes(loadedNodes);

    // Extract LEARNED_PROFILES before replaceState modifies state
    juce::ValueTree profilesTree =
        ProfileStore::findProfilesTree(state, *xmlState);

    if (state.isValid()) {
      // Mark the APVTS swap so parameterChanged can tell a structural change
      // from a session restore apart from live automation (GH #221).
      status.replacingStateForRestore.store(true);
      parameters.replaceState(state);
      status.replacingStateForRestore.store(false);
    }

    if (profilesTree.isValid()) {
      // Serialize the whole restore sequence against profile readers
      // (visualization try-lock in processBlock, getStateInformation)
      juce::SpinLock::ScopedLockType profileLockGuard(profiles.lock());
      profiles.restoreFromTree(profilesTree, engine.group());
    }
  }
}

juce::AudioProcessorEditor* NoiseRepellentAudioProcessor::createEditor() {
  return std::make_unique<NoiseRepellentAudioProcessorEditor>(*this).release();
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
  return std::make_unique<NoiseRepellentAudioProcessor>().release();
}
