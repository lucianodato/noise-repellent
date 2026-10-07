/*
noise-repellent-live -- Zero-Latency Noise Reduction JUCE Plugin

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

NoiseRepellentLiveAudioProcessor::NoiseRepellentLiveAudioProcessor()
    : AudioProcessor(
          BusesProperties()
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      parameters(*this, nullptr, "PARAMETERS", createParameterLayout()),
      dryWetMixer(16384) {
  bypassParameter = dynamic_cast<juce::AudioParameterBool*>(
      parameters.getParameter("bypass"));

  // Route DSP parameter changes through the message thread. bypass is
  // read atomically in processBlock and needs no listener.
  for (const auto& id : {"reduction_amount", "attack_ms", "release_ms",
                         "threshold_db", "knee_db", "learning"}) {
    parameters.addParameterListener(id, this);
  }
}

NoiseRepellentLiveAudioProcessor::~NoiseRepellentLiveAudioProcessor() {
  for (const auto& id : {"reduction_amount", "attack_ms", "release_ms",
                         "threshold_db", "knee_db", "learning"}) {
    parameters.removeParameterListener(id, this);
  }
  cancelPendingUpdate();
  for (auto& engine : engines) {
    engine.reset();
  }
}

juce::AudioProcessorParameter*
NoiseRepellentLiveAudioProcessor::getBypassParameter() const {
  return bypassParameter;
}

juce::AudioProcessorValueTreeState::ParameterLayout
NoiseRepellentLiveAudioProcessor::createParameterLayout() {
  std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "reduction_amount", "Reduction",
      juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 12.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "attack_ms", "Attack", juce::NormalisableRange<float>(0.1f, 100.0f, 0.1f),
      5.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "release_ms", "Release",
      juce::NormalisableRange<float>(10.0f, 1000.0f, 1.0f), 100.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "threshold_db", "Threshold",
      juce::NormalisableRange<float>(-12.0f, 12.0f, 0.1f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "knee_db", "Knee", juce::NormalisableRange<float>(0.0f, 12.0f, 0.1f),
      0.0f));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "bypass", "Internal Bypass", false));

  // Learn is automatable (mirrors RX Voice De-noise "Adaptive Mode") so
  // hosts can drive learn-then-freeze: 1.0 adapts (+ resets the floor on
  // the 0->1 edge), 0.0 freezes. Appended last so existing indices 0-5
  // are unchanged.
  params.push_back(
      std::make_unique<juce::AudioParameterBool>("learning", "Learn", false));

  return {params.begin(), params.end()};
}

void NoiseRepellentLiveAudioProcessor::parameterChanged(
    const juce::String& /*parameterID*/, float /*newValue*/) {
  // Stage the request; the actual engine parameter load happens on the
  // message thread via handleAsyncUpdate (load_parameters is setup-only).
  parametersDirty.store(true, std::memory_order_release);
  triggerAsyncUpdate();
}

void NoiseRepellentLiveAudioProcessor::handleAsyncUpdate() {
  if (parametersDirty.exchange(false, std::memory_order_acq_rel)) {
    applyParameters();
  }
}

void NoiseRepellentLiveAudioProcessor::applyParameters() {
  const juce::ScopedLock sl(getCallbackLock());

  SpecbleachLiveDenoiserParameters p{};
  p.reduction_gain = juce::Decibels::decibelsToGain(
      -parameters.getRawParameterValue("reduction_amount")->load());
  p.attack_time =
      parameters.getRawParameterValue("attack_ms")->load() / 1000.0f;
  p.release_time =
      parameters.getRawParameterValue("release_ms")->load() / 1000.0f;
  // Learn toggle as an automatable parameter: the tracker runs while
  // engaged, then the captured threshold freezes on disengage. The 0->1
  // edge drops the old floor so the tracker converges on the looped
  // noise only. Fast preset so a short noise loop converges.
  const bool wantLearn =
      parameters.getRawParameterValue("learning")->load() > 0.5f;
  if (wantLearn && !learning.load(std::memory_order_acquire)) {
    // Fresh capture (RT-safe flag, honored by process).
    for (auto& engine : engines) {
      if (engine != nullptr) {
        specbleach_live_denoiser_reset_noise_floor(engine.get());
      }
    }
  }
  learning.store(wantLearn, std::memory_order_release);
  p.adaptive_noise = wantLearn;
  p.noise_estimation_method = SPECBLEACH_LIVE_SPP_MMSE;
  p.threshold_db = parameters.getRawParameterValue("threshold_db")->load();
  p.knee_db = parameters.getRawParameterValue("knee_db")->load();

  for (auto& engine : engines) {
    if (engine != nullptr) {
      specbleach_live_denoiser_load_parameters(engine.get(), &p, sizeof(p));
    }
  }
}

void NoiseRepellentLiveAudioProcessor::ensureEnginesInitialized(
    double sampleRate) {
  const uint32_t channels = static_cast<uint32_t>(
      std::max({getTotalNumInputChannels(), getTotalNumOutputChannels(), 1}));
  const bool needRebuild = engines[0] == nullptr ||
                           std::abs(currentSampleRate - sampleRate) > 0.001 ||
                           channels != preparedNumChannels;
  if (!needRebuild) {
    return;
  }

  {
    const juce::ScopedLock sl(getCallbackLock());
    for (auto& engine : engines) {
      engine.reset();
    }
  }

  for (uint32_t ch = 0; ch < channels; ++ch) {
    engines[ch].reset(
        specbleach_live_denoiser_initialize(static_cast<uint32_t>(sampleRate)));
    specbleach_live_denoiser_set_delta_monitoring(
        engines[ch].get(), deltaMonitoring.load(std::memory_order_acquire));
  }

  currentSampleRate = sampleRate;
  preparedNumChannels = channels;
}

void NoiseRepellentLiveAudioProcessor::prepareToPlay(double sampleRate,
                                                     int samplesPerBlock) {
  ensureEnginesInitialized(sampleRate);

  // Zero latency, always. Never call setLatencySamples anywhere else.
  setLatencySamples(0);

  juce::dsp::ProcessSpec spec;
  spec.sampleRate = sampleRate;
  spec.maximumBlockSize =
      static_cast<juce::uint32>(std::max(samplesPerBlock, 16384));
  // DryWetMixer must match the actual bus channel count exactly (its internal
  // delay line asserts on channel mismatch). Hosts re-prepare on layout change.
  spec.numChannels = static_cast<juce::uint32>(
      std::max(getTotalNumInputChannels(), getTotalNumOutputChannels()));
  dryWetMixer.prepare(spec);
  dryWetMixer.setMixingRule(juce::dsp::DryWetMixingRule::linear);
  dryWetMixer.setWetLatency(0.0f);

  applyParameters();

  cancelPendingUpdate();
}

void NoiseRepellentLiveAudioProcessor::releaseResources() {
  dryWetMixer.reset();
}

bool NoiseRepellentLiveAudioProcessor::isBusesLayoutSupported(
    const BusesLayout& layouts) const {
  if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono() &&
      layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo()) {
    return false;
  }

  return layouts.getMainOutputChannelSet() == layouts.getMainInputChannelSet();
}

void NoiseRepellentLiveAudioProcessor::processBlock(
    juce::AudioBuffer<float>& buffer, juce::MidiBuffer&) {
  juce::ScopedNoDenormals noDenormals;
  const int numSamples = buffer.getNumSamples();
  const int numChannels = buffer.getNumChannels();

  if (numSamples == 0 || numChannels == 0) {
    return;
  }

  const bool isBypassed =
      parameters.getRawParameterValue("bypass")->load() > 0.5f;
  const bool monitorDelta = deltaMonitoring.load(std::memory_order_acquire);

  const int procChannels =
      std::min(numChannels, static_cast<int>(engines.size()));

  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);

  // Soft bypass: the engine always processes so gate state and the
  // display stay live while bypassed; the DryWetMixer (~50 ms ramp)
  // crossfades between dry and denoised without clicks.
  for (int ch = 0; ch < procChannels; ++ch) {
    auto& engine = engines[static_cast<size_t>(ch)];
    if (engine != nullptr) {
      specbleach_live_denoiser_process(
          engine.get(), static_cast<uint32_t>(numSamples),
          buffer.getReadPointer(ch), buffer.getWritePointer(ch));
    }
  }

  if (monitorDelta) {
    dryWetMixer.setWetMixProportion(1.0f);
  } else {
    dryWetMixer.setWetMixProportion(isBypassed ? 0.0f : 1.0f);
  }
  dryWetMixer.mixWetSamples(audioBlock);

  // Skip visualization during offline rendering or when the GUI is closed
  if (isNonRealtime() || getActiveEditor() == nullptr) {
    return;
  }

  // Snapshot the engine filterbank (input / post-gate / threshold levels,
  // averaged across channels) for the RX-style band display.
  int start1, size1, start2, size2;
  spectralFifo.prepareToWrite(1, start1, size1, start2, size2);
  if (size1 > 0) {
    BandFrame& frame = spectralBuffer[static_cast<size_t>(start1)];
    frame.inputLevels.fill(0.0f);
    frame.outputLevels.fill(0.0f);
    frame.thresholdLevels.fill(0.0f);

    std::array<float, kNumBands> chIn{};
    std::array<float, kNumBands> chOut{};
    std::array<float, kNumBands> chThr{};
    int numEngines = 0;
    for (int ch = 0; ch < procChannels; ++ch) {
      auto& engine = engines[static_cast<size_t>(ch)];
      if (engine == nullptr) {
        continue;
      }
      if (specbleach_live_denoiser_get_band_levels(
              engine.get(), chIn.data(), chOut.data(), chThr.data())) {
        for (size_t b = 0; b < kNumBands; ++b) {
          frame.inputLevels[b] += chIn[b];
          frame.outputLevels[b] += chOut[b];
          frame.thresholdLevels[b] += chThr[b];
        }
        ++numEngines;
      }
    }
    if (numEngines > 1) {
      const float inv = 1.0f / static_cast<float>(numEngines);
      for (size_t b = 0; b < kNumBands; ++b) {
        frame.inputLevels[b] *= inv;
        frame.outputLevels[b] *= inv;
        frame.thresholdLevels[b] *= inv;
      }
    }
    spectralFifo.finishedWrite(1);
  }
}

void NoiseRepellentLiveAudioProcessor::processBlockBypassed(
    juce::AudioBuffer<float>& buffer, juce::MidiBuffer&) {
  const int numSamples = buffer.getNumSamples();
  const int numChannels = buffer.getNumChannels();

  if (numSamples == 0 || numChannels == 0) {
    return;
  }

  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);
  dryWetMixer.setWetMixProportion(0.0f);
  dryWetMixer.mixWetSamples(audioBlock);
}

bool NoiseRepellentLiveAudioProcessor::getNextBandFrame(BandFrame& frame) {
  int start1, size1, start2, size2;
  spectralFifo.prepareToRead(1, start1, size1, start2, size2);

  if (size1 > 0) {
    frame = spectralBuffer[static_cast<size_t>(start1)];
    spectralFifo.finishedRead(1);
    return true;
  }
  return false;
}

void NoiseRepellentLiveAudioProcessor::setLearning(bool shouldLearn) {
  // Route through the automatable parameter so hosts, state, and the GUI
  // button stay in sync; the floor reset happens on the 0->1 edge in
  // applyParameters().
  if (auto* param = parameters.getParameter("learning")) {
    param->setValueNotifyingHost(shouldLearn ? 1.0f : 0.0f);
  }
}

void NoiseRepellentLiveAudioProcessor::setDeltaMonitoring(bool shouldMonitor) {
  deltaMonitoring.store(shouldMonitor, std::memory_order_release);
  for (auto& engine : engines) {
    if (engine != nullptr) {
      specbleach_live_denoiser_set_delta_monitoring(engine.get(),
                                                    shouldMonitor);
    }
  }
}

bool NoiseRepellentLiveAudioProcessor::getLiveBandEdges(float* lowerHz,
                                                        float* upperHz) {
  for (auto& engine : engines) {
    if (engine != nullptr) {
      return specbleach_live_denoiser_get_band_edges(engine.get(), lowerHz,
                                                     upperHz);
    }
  }
  return false;
}

juce::AudioProcessorEditor* NoiseRepellentLiveAudioProcessor::createEditor() {
  return new NoiseRepellentLiveAudioProcessorEditor(*this);
}

void NoiseRepellentLiveAudioProcessor::getStateInformation(
    juce::MemoryBlock& destData) {
  auto state = parameters.copyState();
  std::unique_ptr<juce::XmlElement> xml(state.createXml());
  copyXmlToBinary(*xml, destData);
}

void NoiseRepellentLiveAudioProcessor::setStateInformation(const void* data,
                                                           int sizeInBytes) {
  std::unique_ptr<juce::XmlElement> xmlState(
      getXmlFromBinary(data, sizeInBytes));
  if (xmlState != nullptr) {
    parameters.replaceState(juce::ValueTree::fromXml(*xmlState));
  }
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
  return new NoiseRepellentLiveAudioProcessor();
}
