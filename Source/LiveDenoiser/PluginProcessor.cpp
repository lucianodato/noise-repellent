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
#include <cstring>

NoiseRepellentLiveAudioProcessor::NoiseRepellentLiveAudioProcessor()
    : AudioProcessor(
          BusesProperties()
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      parameters(*this, nullptr, "PARAMETERS", createParameterLayout()),
      dryWetMixer(16384) {
  bypassParameter = dynamic_cast<juce::AudioParameterBool*>(
      parameters.getParameter("bypass"));
  adaptiveMethodParameter = dynamic_cast<juce::AudioParameterChoice*>(
      parameters.getParameter("adaptive_method"));

  // Route DSP parameter changes through the message thread. bypass is
  // read atomically in processBlock and needs no listener.
  for (const auto& id :
       {"reduction_amount", "adaptive_noise", "adaptive_method",
        "smoothing_factor", "suppression_strength"}) {
    parameters.addParameterListener(id, this);
  }
}

NoiseRepellentLiveAudioProcessor::~NoiseRepellentLiveAudioProcessor() {
  for (const auto& id :
       {"reduction_amount", "adaptive_noise", "adaptive_method",
        "smoothing_factor", "suppression_strength"}) {
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

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "adaptive_noise", "Adaptive Noise", true));

  params.push_back(std::make_unique<juce::AudioParameterChoice>(
      "adaptive_method", "Estimation Method",
      juce::StringArray{"SPP-MMSE (Unbiased)", "Brandt (Trimmed Mean)",
                        "Martin (Min Statistics)"},
      2));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "smoothing_factor", "Smoothing",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 0.0f));

  params.push_back(std::make_unique<juce::AudioParameterFloat>(
      "suppression_strength", "Aggressiveness",
      juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 50.0f));

  params.push_back(std::make_unique<juce::AudioParameterBool>(
      "bypass", "Internal Bypass", false));

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

  SpecbleachTimeDenoiserParameters p{};
  p.reduction_gain = juce::Decibels::decibelsToGain(
      -parameters.getRawParameterValue("reduction_amount")->load());
  p.smoothing_factor =
      parameters.getRawParameterValue("smoothing_factor")->load() / 100.0f;
  p.adaptive_noise =
      parameters.getRawParameterValue("adaptive_noise")->load() > 0.5f;
  p.noise_estimation_method =
      static_cast<SpecbleachNoiseEstimationMethod>(
          adaptiveMethodParameter != nullptr
              ? adaptiveMethodParameter->getIndex()
              : 2);
  p.suppression_strength =
      parameters.getRawParameterValue("suppression_strength")->load() / 100.0f;

  for (auto& engine : engines) {
    if (engine != nullptr) {
      specbleach_time_denoiser_load_parameters(engine.get(), &p, sizeof(p));
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
        specbleach_time_denoiser_initialize(static_cast<uint32_t>(sampleRate)));
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

  // Pre-allocate persistent buffers to prevent audio-thread allocations
  dryInputL.resize(static_cast<size_t>(std::max(samplesPerBlock, 16384)), 0.0f);

  applyParameters();

  // Reset FFT accumulation
  fftAccumInput.fill(0.0f);
  fftAccumOutput.fill(0.0f);
  fftAccumCount = 0;
  silenceVisualCounter = 0;
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

  // Save dry input copy for FFT visualization before in-place processing
  const size_t copySamples =
      std::min(static_cast<size_t>(numSamples), dryInputL.size());
  if (numChannels >= 1 && copySamples > 0) {
    std::copy_n(buffer.getReadPointer(0), copySamples, dryInputL.begin());
  }

  juce::dsp::AudioBlock<float> audioBlock(buffer);
  dryWetMixer.pushDrySamples(audioBlock);

  if (!isBypassed) {
    const int procChannels =
        std::min(numChannels, static_cast<int>(engines.size()));
    for (int ch = 0; ch < procChannels; ++ch) {
      auto& engine = engines[static_cast<size_t>(ch)];
      if (engine != nullptr) {
        specbleach_time_denoiser_process(
            engine.get(), static_cast<uint32_t>(numSamples),
            buffer.getReadPointer(ch), buffer.getWritePointer(ch));
      }
    }
  }

  // Soft crossfade bypass using JUCE DryWetMixer (wet latency is 0)
  dryWetMixer.setWetMixProportion(isBypassed ? 0.0f : 1.0f);
  dryWetMixer.mixWetSamples(audioBlock);

  // Skip FFT analysis during offline rendering or when the GUI is closed
  if (isNonRealtime() || getActiveEditor() == nullptr) {
    fftAccumCount = 0;
    return;
  }

  // Detect silence for cheap visualization frames
  bool isSilent = true;
  for (int ch = 0; ch < numChannels && isSilent; ++ch) {
    const float* d = buffer.getReadPointer(ch);
    for (int i = 0; i < numSamples; ++i) {
      if (std::abs(d[i]) > 1e-5f) {
        isSilent = false;
        break;
      }
    }
  }

  if (isSilent) {
    fftAccumCount = 0;
    if ((++silenceVisualCounter % 4) == 0) {
      int start1, size1, start2, size2;
      spectralFifo.prepareToWrite(1, start1, size1, start2, size2);
      if (size1 > 0) {
        SpectralFrame& frame = spectralBuffer[static_cast<size_t>(start1)];
        frame.inputMagnitudeDB.fill(-120.0f);
        frame.outputMagnitudeDB.fill(-120.0f);
        spectralFifo.finishedWrite(1);
      }
    }
    return;
  }

  // Accumulate samples until a full FFT window (hop = 75%)
  const float* inputSrc = dryInputL.data();
  const float* outputSrc = buffer.getReadPointer(0);
  for (int s = 0; s < numSamples && static_cast<size_t>(s) < copySamples; ++s) {
    if (fftAccumCount < kFftSize) {
      fftAccumInput[fftAccumCount] = inputSrc[s];
      fftAccumOutput[fftAccumCount] = outputSrc[s];
      fftAccumCount++;
    }

    if (fftAccumCount >= kFftSize) {
      int start1, size1, start2, size2;
      spectralFifo.prepareToWrite(1, start1, size1, start2, size2);

      if (size1 > 0) {
        SpectralFrame& frame = spectralBuffer[static_cast<size_t>(start1)];

        std::memcpy(fftInputWork.data(), fftAccumInput.data(),
                    kFftSize * sizeof(float));
        std::fill(fftInputWork.begin() + kFftSize, fftInputWork.end(), 0.0f);
        fftWindow.multiplyWithWindowingTable(fftInputWork.data(), kFftSize);
        fftAnalyzer.performFrequencyOnlyForwardTransform(fftInputWork.data());
        for (size_t i = 0; i < kFftBins; ++i) {
          const float mag = fftInputWork[i] / static_cast<float>(kFftBins);
          frame.inputMagnitudeDB[i] = 20.0f * std::log10(std::max(mag, 1e-7f));
        }

        std::memcpy(fftOutputWork.data(), fftAccumOutput.data(),
                    kFftSize * sizeof(float));
        std::fill(fftOutputWork.begin() + kFftSize, fftOutputWork.end(), 0.0f);
        fftWindow.multiplyWithWindowingTable(fftOutputWork.data(), kFftSize);
        fftAnalyzer.performFrequencyOnlyForwardTransform(fftOutputWork.data());
        for (size_t i = 0; i < kFftBins; ++i) {
          const float mag = fftOutputWork[i] / static_cast<float>(kFftBins);
          frame.outputMagnitudeDB[i] = 20.0f * std::log10(std::max(mag, 1e-7f));
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

bool NoiseRepellentLiveAudioProcessor::getNextSpectralFrame(
    SpectralFrame& frame) {
  int start1, size1, start2, size2;
  spectralFifo.prepareToRead(1, start1, size1, start2, size2);

  if (size1 > 0) {
    frame = spectralBuffer[static_cast<size_t>(start1)];
    spectralFifo.finishedRead(1);
    return true;
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
