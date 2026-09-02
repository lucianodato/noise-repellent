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

#include <cmath>
#include <random>

#include <juce_core/juce_core.h>
#include <juce_events/juce_events.h>
#include <juce_gui_basics/juce_gui_basics.h>

#include "../LiveDenoiser/PluginProcessor.h"

namespace {

void pumpMessageLoop(int maxMillis = 50) {
  if (auto* mm = juce::MessageManager::getInstanceWithoutCreating()) {
    mm->runDispatchLoopUntil(maxMillis);
  }
}

void generateNoiseBuffer(juce::AudioBuffer<float>& buffer, float rms = 0.05f) {
  std::mt19937 gen(1337);
  std::normal_distribution<float> dist(0.0f, rms);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] = dist(gen);
    }
  }
}

void setParam(NoiseRepellentLiveAudioProcessor& proc,
              const juce::String& paramId, float value) {
  if (auto* param = proc.getAPVTS().getParameter(paramId)) {
    param->setValueNotifyingHost(param->convertTo0to1(value));
  }
  pumpMessageLoop(5);
}

float bufferRms(const juce::AudioBuffer<float>& buffer) {
  double sum = 0.0;
  size_t count = 0;
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    const auto* data = buffer.getReadPointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      sum += static_cast<double>(data[s]) * static_cast<double>(data[s]);
      ++count;
    }
  }
  return count > 0
             ? static_cast<float>(std::sqrt(sum / static_cast<double>(count)))
             : 0.0f;
}

} // namespace

class LiveUxTest : public juce::UnitTest {
public:
  LiveUxTest() : juce::UnitTest("Live Plugin UX Test", "NoiseRepellentLive") {
  }

  void runTest() override {
    beginTest("zero latency invariant");
    {
      NoiseRepellentLiveAudioProcessor proc;
      expectEquals(proc.getLatencySamples(), 0);

      proc.setRateAndBufferSizeDetails(48000.0, 512);
      proc.prepareToPlay(48000.0, 512);
      expectEquals(proc.getLatencySamples(), 0);

      juce::AudioBuffer<float> buffer(2, 512);
      juce::MidiBuffer midi;
      generateNoiseBuffer(buffer);
      proc.processBlock(buffer, midi);
      expectEquals(proc.getLatencySamples(), 0);

      proc.releaseResources();
      expectEquals(proc.getLatencySamples(), 0);
    }

    beginTest("parameter set and defaults");
    {
      NoiseRepellentLiveAudioProcessor proc;
      auto& apvts = proc.getAPVTS();

      const char* expectedIds[] = {"reduction_amount",     "adaptive_noise",
                                   "adaptive_method",      "smoothing_factor",
                                   "suppression_strength", "bypass"};
      for (const auto* id : expectedIds) {
        expect(apvts.getParameter(id) != nullptr,
               juce::String("missing parameter: ") + id);
      }

      // No stale parameters from the main plugin
      const int numParams = static_cast<int>(proc.getParameters().size());
      expectEquals(numParams, 6);

      auto* reduction = static_cast<juce::AudioParameterFloat*>(
          apvts.getParameter("reduction_amount"));
      expectWithinAbsoluteError(reduction->get(), 12.0f, 1e-5f);
      expectWithinAbsoluteError(reduction->range.start, 0.0f, 1e-5f);
      expectWithinAbsoluteError(reduction->range.end, 40.0f, 1e-5f);

      auto* adaptive = static_cast<juce::AudioParameterBool*>(
          apvts.getParameter("adaptive_noise"));
      expect(adaptive->get());

      auto* bypass =
          static_cast<juce::AudioParameterBool*>(apvts.getParameter("bypass"));
      expect(!bypass->get());
    }

    beginTest("bus layout support");
    {
      NoiseRepellentLiveAudioProcessor proc;
      juce::AudioProcessor::BusesLayout mono;
      mono.inputBuses.add(juce::AudioChannelSet::mono());
      mono.outputBuses.add(juce::AudioChannelSet::mono());
      expect(proc.isBusesLayoutSupported(mono));

      juce::AudioProcessor::BusesLayout stereo;
      stereo.inputBuses.add(juce::AudioChannelSet::stereo());
      stereo.outputBuses.add(juce::AudioChannelSet::stereo());
      expect(proc.isBusesLayoutSupported(stereo));

      juce::AudioProcessor::BusesLayout mismatched;
      mismatched.inputBuses.add(juce::AudioChannelSet::mono());
      mismatched.outputBuses.add(juce::AudioChannelSet::stereo());
      expect(!proc.isBusesLayoutSupported(mismatched));

      juce::AudioProcessor::BusesLayout surround;
      surround.inputBuses.add(juce::AudioChannelSet::quadraphonic());
      surround.outputBuses.add(juce::AudioChannelSet::quadraphonic());
      expect(!proc.isBusesLayoutSupported(surround));
    }

    beginTest("real-time processing across odd block sizes");
    {
      NoiseRepellentLiveAudioProcessor proc;
      proc.setRateAndBufferSizeDetails(48000.0, 512);
      proc.prepareToPlay(48000.0, 512);

      juce::MidiBuffer midi;
      const int blockSizes[] = {1, 13, 64, 512};
      for (int blockSize : blockSizes) {
        juce::AudioBuffer<float> buffer(2, blockSize);
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
        for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
          const auto* data = buffer.getReadPointer(ch);
          for (int s = 0; s < blockSize; ++s) {
            expect(std::isfinite(data[s]), "non-finite sample in output");
          }
        }
      }
    }

    beginTest("reduction reduces noise rms");
    {
      // Regression: block sizes above the engine's internal ring capacity
      // used to starve analysis entirely (silent pass-through)
      for (const int blockSize : {512, 4096}) {
        NoiseRepellentLiveAudioProcessor proc;
        proc.setRateAndBufferSizeDetails(48000.0, blockSize);
        proc.prepareToPlay(48000.0, blockSize);

        // Warm up in adaptive mode so the engine learns the noise floor
        juce::MidiBuffer midi;
        juce::AudioBuffer<float> buffer(2, blockSize);
        for (int i = 0; i < 96; ++i) {
          generateNoiseBuffer(buffer);
          proc.processBlock(buffer, midi);
        }
        pumpMessageLoop(10);

        // Freeze the learned estimate and max out reduction
        setParam(proc, "adaptive_noise", 0.0f);
        setParam(proc, "reduction_amount", 40.0f);
        setParam(proc, "suppression_strength", 100.0f);
        for (int i = 0; i < 4; ++i) {
          generateNoiseBuffer(buffer);
          proc.processBlock(buffer, midi);
        }

        // Measure dry RMS of the same noise distribution
        juce::AudioBuffer<float> dry(2, 512 * 16);
        generateNoiseBuffer(dry);
        const float dryRms = bufferRms(dry);

        // Measure wet RMS
        juce::AudioBuffer<float> wet(2, 512 * 16);
        for (int offset = 0; offset < wet.getNumSamples(); offset += 512) {
          generateNoiseBuffer(buffer);
          proc.processBlock(buffer, midi);
          for (int ch = 0; ch < 2; ++ch) {
            wet.copyFrom(ch, offset, buffer, ch, 0, 512);
          }
        }
        const float wetRms = bufferRms(wet);

        expect(wetRms < dryRms * 0.3f,
               "expected strong noise reduction at max settings, block " +
                   juce::String(blockSize));

        proc.releaseResources();
      }
    }

    beginTest("state round trip");
    {
      NoiseRepellentLiveAudioProcessor proc;
      setParam(proc, "reduction_amount", 25.0f);
      setParam(proc, "adaptive_noise", 0.0f);
      setParam(proc, "smoothing_factor", 60.0f);
      setParam(proc, "suppression_strength", 80.0f);
      setParam(proc, "adaptive_method", 0.0f);

      juce::MemoryBlock state;
      proc.getStateInformation(state);

      NoiseRepellentLiveAudioProcessor proc2;
      proc2.setStateInformation(state.getData(),
                                static_cast<int>(state.getSize()));

      auto& apvts2 = proc2.getAPVTS();
      expectWithinAbsoluteError(static_cast<juce::AudioParameterFloat*>(
                                    apvts2.getParameter("reduction_amount"))
                                    ->get(),
                                25.0f, 1e-5f);
      expect(!static_cast<juce::AudioParameterBool*>(
                  apvts2.getParameter("adaptive_noise"))
                  ->get());
      expectWithinAbsoluteError(static_cast<juce::AudioParameterFloat*>(
                                    apvts2.getParameter("smoothing_factor"))
                                    ->get(),
                                60.0f, 1e-5f);
      expectWithinAbsoluteError(static_cast<juce::AudioParameterFloat*>(
                                    apvts2.getParameter("suppression_strength"))
                                    ->get(),
                                80.0f, 1e-5f);
      expectEquals(static_cast<juce::AudioParameterChoice*>(
                       apvts2.getParameter("adaptive_method"))
                       ->getIndex(),
                   0);
    }
  }
};

static LiveUxTest liveUxTest;

int main(int argc, char* argv[]) {
  juce::ScopedJuceInitialiser_GUI guiInit;
  juce::ignoreUnused(argc, argv);

  juce::UnitTestRunner runner;
  runner.setAssertOnFailure(false);
  runner.runAllTests();

  int numFailed = 0;
  for (int i = 0; i < runner.getNumResults(); ++i) {
    if (const auto* r = runner.getResult(i)) {
      numFailed += r->failures;
    }
  }
  return numFailed > 0 ? 1 : 0;
}
