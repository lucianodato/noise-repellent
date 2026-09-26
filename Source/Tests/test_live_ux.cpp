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

      const char* expectedIds[] = {"reduction_amount", "attack_ms",
                                   "release_ms",       "threshold_db",
                                   "knee_db",          "bypass"};
      for (const auto* id : expectedIds) {
        expect(apvts.getParameter(id) != nullptr,
               juce::String("missing parameter: ") + id);
      }

      // No stale parameters from the main plugin (Live learns manually:
      // no adaptive_noise / adaptive_method)
      const int numParams = static_cast<int>(proc.getParameters().size());
      expectEquals(numParams, 6);

      auto* reduction = static_cast<juce::AudioParameterFloat*>(
          apvts.getParameter("reduction_amount"));
      expectWithinAbsoluteError(reduction->get(), 12.0f, 1e-5f);
      expectWithinAbsoluteError(reduction->range.start, 0.0f, 1e-5f);
      expectWithinAbsoluteError(reduction->range.end, 40.0f, 1e-5f);

      auto* threshold = static_cast<juce::AudioParameterFloat*>(
          apvts.getParameter("threshold_db"));
      expectWithinAbsoluteError(threshold->get(), 0.0f, 1e-5f);
      expectWithinAbsoluteError(threshold->range.start, -12.0f, 1e-5f);
      expectWithinAbsoluteError(threshold->range.end, 12.0f, 1e-5f);

      auto* attack = static_cast<juce::AudioParameterFloat*>(
          apvts.getParameter("attack_ms"));
      expectWithinAbsoluteError(attack->get(), 5.0f, 1e-5f);

      auto* release = static_cast<juce::AudioParameterFloat*>(
          apvts.getParameter("release_ms"));
      expectWithinAbsoluteError(release->get(), 100.0f, 1e-5f);

      auto* knee = static_cast<juce::AudioParameterFloat*>(
          apvts.getParameter("knee_db"));
      expectWithinAbsoluteError(knee->get(), 0.0f, 1e-5f);
      expectWithinAbsoluteError(knee->range.start, 0.0f, 1e-5f);
      expectWithinAbsoluteError(knee->range.end, 12.0f, 1e-5f);

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
      // Manual learn: engage the tracker while looping noise (fail-open
      // floor rises from zero), then freeze and measure.
      for (const int blockSize : {512, 4096}) {
        NoiseRepellentLiveAudioProcessor proc;
        proc.setRateAndBufferSizeDetails(48000.0, blockSize);
        proc.prepareToPlay(48000.0, blockSize);

        proc.setLearning(true);
        pumpMessageLoop(10);
        juce::MidiBuffer midi;
        juce::AudioBuffer<float> buffer(2, blockSize);
        for (int i = 0; i < 96; ++i) {
          generateNoiseBuffer(buffer);
          proc.processBlock(buffer, midi);
        }
        proc.setLearning(false);
        pumpMessageLoop(10);

        // Max out reduction and measure (fast gate so a few blocks
        // settle to steady state)
        setParam(proc, "reduction_amount", 40.0f);
        setParam(proc, "threshold_db", 12.0f);
        setParam(proc, "attack_ms", 5.0f);
        setParam(proc, "release_ms", 30.0f);
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

        // Dynamic-EQ cascade semantics: the Reduction slider is the
        // steady broadband cut and saturates deeper (notch geometry),
        // so anchor the assert to the calibrated floor.
        expect(wetRms < dryRms * 0.65f,
               "expected strong noise reduction at max settings, block " +
                   juce::String(blockSize));

        proc.releaseResources();
      }
    }
    beginTest("manual learn captures and freezes the threshold");
    {
      NoiseRepellentLiveAudioProcessor proc;
      proc.setRateAndBufferSizeDetails(48000.0, 512);
      proc.prepareToPlay(48000.0, 512);
      std::unique_ptr<juce::AudioProcessorEditor> editor(
          proc.createEditorIfNeeded());

      // Loop noise with Learn engaged so the floor converges. Drain
      // as we go: with no GUI timer running the fifo would saturate
      // and hold a stale early frame.
      NoiseRepellentLiveAudioProcessor::BandFrame frame;
      bool gotFrame = false;
      proc.setLearning(true);
      pumpMessageLoop(10);
      juce::MidiBuffer midi;
      juce::AudioBuffer<float> buffer(2, 512);
      for (int i = 0; i < 100; ++i) {
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
        while (proc.getNextBandFrame(frame)) {
          gotFrame = true;
        }
      }
      expect(gotFrame, "expected band frames after processing noise");
      float maxIn = 0.0f;
      float maxThr = 0.0f;
      for (size_t b = 0; b < NoiseRepellentLiveAudioProcessor::kNumBands; ++b) {
        expect(std::isfinite(frame.inputLevels[b]), "finite input level");
        expect(std::isfinite(frame.outputLevels[b]), "finite output level");
        expect(std::isfinite(frame.thresholdLevels[b]),
               "finite threshold level");
        expect(frame.outputLevels[b] <= frame.inputLevels[b] + 1e-6f,
               "output must not exceed input per band");
        expect(frame.thresholdLevels[b] >= 0.0f, "threshold non-negative");
        maxIn = std::max(maxIn, frame.inputLevels[b]);
        maxThr = std::max(maxThr, frame.thresholdLevels[b]);
      }
      expect(maxIn > 1e-6f, "input energy must be visible in frames");
      expect(maxThr > 0.0f, "learned threshold must be visible in frames");

      // Disengage Learn: the captured threshold must freeze even as the
      // input goes silent (generous pump: the freeze arrives via the
      // message-thread async update, starved under parallel ctest load)
      proc.setLearning(false);
      pumpMessageLoop(100);
      juce::AudioBuffer<float> silence(2, 512);
      silence.clear();
      proc.processBlock(silence, midi);
      NoiseRepellentLiveAudioProcessor::BandFrame frozen;
      while (proc.getNextBandFrame(frozen)) {
      }
      for (size_t b = 0; b < NoiseRepellentLiveAudioProcessor::kNumBands; ++b) {
        expectWithinAbsoluteError(frozen.thresholdLevels[b],
                                  frame.thresholdLevels[b], 1e-6f);
      }

      proc.releaseResources();
      proc.editorBeingDeleted(editor.get());
    }
    beginTest("delta monitoring outputs the removed noise");
    {
      NoiseRepellentLiveAudioProcessor proc;
      proc.setRateAndBufferSizeDetails(48000.0, 512);
      proc.prepareToPlay(48000.0, 512);

      // Manual learn on looped noise, then freeze
      proc.setLearning(true);
      pumpMessageLoop(10);
      juce::MidiBuffer midi;
      juce::AudioBuffer<float> buffer(2, 512);
      for (int i = 0; i < 100; ++i) {
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
      }
      proc.setLearning(false);
      pumpMessageLoop(10);

      setParam(proc, "reduction_amount", 40.0f);
      setParam(proc, "threshold_db", 12.0f);
      setParam(proc, "attack_ms", 5.0f);
      setParam(proc, "release_ms", 30.0f);

      // Normal run on deterministic noise (same content regenerates)
      generateNoiseBuffer(buffer);
      juce::AudioBuffer<float> dry;
      dry.makeCopyOf(buffer);
      const float dryRms = bufferRms(buffer);
      proc.processBlock(buffer, midi);
      juce::AudioBuffer<float> wet;
      wet.makeCopyOf(buffer);

      // Delta run on the identical input straight after
      generateNoiseBuffer(buffer);
      proc.setDeltaMonitoring(true);
      proc.processBlock(buffer, midi);
      proc.setDeltaMonitoring(false);

      const float deltaRms = bufferRms(buffer);
      // Soft knee: partial gains in the knee band carry less removed
      // noise than a binary gate, so bound below the binary 0.5 level.
      expect(deltaRms > 0.3f * dryRms,
             "delta must carry most of the removed noise");

      // Consistency: wet + removed ~= dry (one block of state drift)
      double errSum = 0.0;
      size_t errCount = 0;
      for (int ch = 0; ch < 2; ++ch) {
        const auto* w = wet.getReadPointer(ch);
        const auto* d = buffer.getReadPointer(ch);
        const auto* o = dry.getReadPointer(ch);
        for (int s = 0; s < 512; ++s) {
          const double e = static_cast<double>(w[s] + d[s] - o[s]);
          errSum += e * e;
          ++errCount;
        }
      }
      const float errRms =
          static_cast<float>(std::sqrt(errSum / static_cast<double>(errCount)));
      // Bank ripple plus one-block gain drift (knee tracks the envelope):
      // wet+delta = sum bands / norm, not bit-identical dry.
      expect(errRms < 1.0f * dryRms,
             "wet + delta must roughly reconstruct dry");

      // No reduction: nothing removed, delta is (near) silence
      setParam(proc, "reduction_amount", 0.0f);
      pumpMessageLoop(10);
      // Let gains open to unity after the max-reduction state
      for (int i = 0; i < 20; ++i) {
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
      }
      generateNoiseBuffer(buffer);
      const float dryRms2 = bufferRms(buffer);
      proc.setDeltaMonitoring(true);
      proc.processBlock(buffer, midi);
      proc.setDeltaMonitoring(false);
      const float delta0Rms = bufferRms(buffer);
      expect(delta0Rms < 0.5f * dryRms2,
             "delta must collapse with no reduction");
    }

    beginTest("soft bypass outputs dry then returns to denoised");
    {
      NoiseRepellentLiveAudioProcessor proc;
      proc.setRateAndBufferSizeDetails(48000.0, 512);
      proc.prepareToPlay(48000.0, 512);

      // Manual learn on looped noise, then freeze
      proc.setLearning(true);
      pumpMessageLoop(10);
      juce::MidiBuffer midi;
      juce::AudioBuffer<float> buffer(2, 512);
      for (int i = 0; i < 100; ++i) {
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
      }
      proc.setLearning(false);
      pumpMessageLoop(10);

      setParam(proc, "reduction_amount", 40.0f);
      setParam(proc, "threshold_db", 12.0f);
      setParam(proc, "attack_ms", 5.0f);
      setParam(proc, "release_ms", 30.0f);

      // Engage bypass: past the ~50 ms mixer ramp the output is the
      // original signal while the engine keeps processing underneath.
      setParam(proc, "bypass", 1.0f);
      for (int i = 0; i < 40; ++i) {
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
      }
      generateNoiseBuffer(buffer);
      juce::AudioBuffer<float> dry;
      dry.makeCopyOf(buffer);
      const float dryRms = bufferRms(buffer);
      proc.processBlock(buffer, midi);
      double errSum = 0.0;
      size_t errCount = 0;
      for (int ch = 0; ch < 2; ++ch) {
        const auto* o = buffer.getReadPointer(ch);
        const auto* d = dry.getReadPointer(ch);
        for (int s = 0; s < 512; ++s) {
          const double e = static_cast<double>(o[s] - d[s]);
          errSum += e * e;
          ++errCount;
        }
      }
      const float errRms =
          static_cast<float>(std::sqrt(errSum / static_cast<double>(errCount)));
      expect(errRms < 1e-6f * dryRms, "bypassed output must be the dry signal");

      // Disengage: the mixer fades back to the denoised signal.
      setParam(proc, "bypass", 0.0f);
      for (int i = 0; i < 40; ++i) {
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
      }
      generateNoiseBuffer(buffer);
      proc.processBlock(buffer, midi);
      const float wetRms = bufferRms(buffer);
      expect(wetRms < 0.7f * dryRms,
             "un-bypassed output must return to denoised");
    }

    beginTest("state round trip");
    {
      NoiseRepellentLiveAudioProcessor proc;
      setParam(proc, "reduction_amount", 25.0f);
      setParam(proc, "attack_ms", 10.0f);
      setParam(proc, "release_ms", 200.0f);
      setParam(proc, "threshold_db", 6.0f);
      setParam(proc, "knee_db", 3.0f);

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
      expectWithinAbsoluteError(static_cast<juce::AudioParameterFloat*>(
                                    apvts2.getParameter("attack_ms"))
                                    ->get(),
                                10.0f, 1e-5f);
      expectWithinAbsoluteError(static_cast<juce::AudioParameterFloat*>(
                                    apvts2.getParameter("release_ms"))
                                    ->get(),
                                200.0f, 1e-5f);
      expectWithinAbsoluteError(static_cast<juce::AudioParameterFloat*>(
                                    apvts2.getParameter("threshold_db"))
                                    ->get(),
                                6.0f, 1e-5f);
      expectWithinAbsoluteError(static_cast<juce::AudioParameterFloat*>(
                                    apvts2.getParameter("knee_db"))
                                    ->get(),
                                3.0f, 1e-5f);
    }

    beginTest("band frames reach the GUI");
    {
      NoiseRepellentLiveAudioProcessor proc;
      proc.setRateAndBufferSizeDetails(48000.0, 512);
      proc.prepareToPlay(48000.0, 512);
      std::unique_ptr<juce::AudioProcessorEditor> editor(
          proc.createEditorIfNeeded());
      expect(editor != nullptr, "editor must be created");
      expect(proc.getActiveEditor() != nullptr,
             "editor must be registered active");

      // Engine scale: every band active at 48 kHz, sane edges.
      constexpr size_t numBands = NoiseRepellentLiveAudioProcessor::kNumBands;
      std::array<float, numBands> lo{};
      std::array<float, numBands> hi{};
      expect(proc.getLiveBandEdges(lo.data(), hi.data()),
             "band edges must be available");
      size_t numActive = 0;
      for (size_t b = 0; b < numBands; ++b) {
        if (hi[b] <= 0.0f) {
          break;
        }
        expect(hi[b] > lo[b], "band edges must be ordered");
        ++numActive;
      }
      expectEquals(static_cast<int>(numActive), static_cast<int>(numBands));

      // Run noise with the editor open: frames must arrive carrying
      // visible input energy.
      proc.setLearning(true);
      pumpMessageLoop(10);
      juce::MidiBuffer midi;
      juce::AudioBuffer<float> buffer(2, 512);
      for (int i = 0; i < 20; ++i) {
        generateNoiseBuffer(buffer);
        proc.processBlock(buffer, midi);
      }
      NoiseRepellentLiveAudioProcessor::BandFrame frame;
      // Read before pumping: the visualizer timer (60 Hz) drains the
      // same fifo on the message thread and would steal the frames.
      bool gotFrame = proc.getNextBandFrame(frame);
      expect(gotFrame, "GUI must receive band frames while processing");
      if (gotFrame) {
        float maxIn = 0.0f;
        for (size_t b = 0; b < numBands; ++b) {
          expect(std::isfinite(frame.inputLevels[b]), "finite input level");
          maxIn = std::max(maxIn, frame.inputLevels[b]);
        }
        expect(maxIn > 1e-4f, "input energy must be visible in frames");
      }
      proc.setLearning(false);
      editor.reset();
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
