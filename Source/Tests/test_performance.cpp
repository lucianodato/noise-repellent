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

#include <chrono>
#include <random>

#include <juce_core/juce_core.h>
#include <juce_events/juce_events.h>
#include <juce_gui_basics/juce_gui_basics.h>

#include "PluginProcessor.h"

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

void setParam(NoiseRepellentAudioProcessor& proc, const juce::String& paramId,
              float value) {
  if (auto* param = proc.getAPVTS().getParameter(paramId)) {
    param->setValueNotifyingHost(param->convertTo0to1(value));
  }
  pumpMessageLoop(5);
}

} // namespace

class PerformanceRatioTest : public juce::UnitTest {
public:
  PerformanceRatioTest()
      : juce::UnitTest("Performance Ratio Tests", "NoiseRepellent") {
  }

  void runTest() override {
    beginTest("Silent Idle vs Active Denoising Throughput Ratio");
    testSilentIdleThroughputRatio();

    beginTest("Bypass vs Active Denoising Throughput Ratio");
    testBypassThroughputRatio();
  }

private:
  void testSilentIdleThroughputRatio() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    constexpr int numBlocks = 600;

    // 1. Measure active denoising duration
    int64_t activeDurationUs = 0;
    {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(20);

      juce::AudioBuffer<float> buffer(2, blockSize);
      juce::MidiBuffer midi;

      // Learn profile
      setParam(proc, "learn_noise", 1.0f);
      for (int i = 0; i < 50; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }
      setParam(proc, "learn_noise", 0.0f);
      pumpMessageLoop(20);

      // Warmup
      for (int i = 0; i < 50; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }

      // Timed active processing
      auto start = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < numBlocks; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }
      auto end = std::chrono::high_resolution_clock::now();
      activeDurationUs =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();

      proc.releaseResources();
    }

    // 2. Measure silent idle duration (no profile, silent input)
    int64_t silentDurationUs = 0;
    {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(20);

      juce::AudioBuffer<float> buffer(2, blockSize);
      juce::MidiBuffer midi;

      // Warmup
      for (int i = 0; i < 50; ++i) {
        buffer.clear();
        proc.processBlock(buffer, midi);
      }

      // Timed silent processing
      auto start = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < numBlocks; ++i) {
        buffer.clear();
        proc.processBlock(buffer, midi);
      }
      auto end = std::chrono::high_resolution_clock::now();
      silentDurationUs =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();

      proc.releaseResources();
    }

    expect(silentDurationUs > 0, "Silent duration must be greater than zero");
    expect(activeDurationUs > 0, "Active duration must be greater than zero");

    const double speedupRatio =
        static_cast<double>(activeDurationUs) /
        static_cast<double>(std::max<int64_t>(1, silentDurationUs));

    expect(speedupRatio >= 2.5,
           "Silent idle processing must be at least 2.5x faster than active "
           "denoising");
  }

  void testBypassThroughputRatio() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    constexpr int numBlocks = 600;

    int64_t activeDurationUs = 0;
    {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(20);

      juce::AudioBuffer<float> buffer(2, blockSize);
      juce::MidiBuffer midi;

      setParam(proc, "learn_noise", 1.0f);
      for (int i = 0; i < 50; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }
      setParam(proc, "learn_noise", 0.0f);
      pumpMessageLoop(20);

      for (int i = 0; i < 50; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }

      auto start = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < numBlocks; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }
      auto end = std::chrono::high_resolution_clock::now();
      activeDurationUs =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();

      proc.releaseResources();
    }

    int64_t bypassDurationUs = 0;
    {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(20);

      juce::AudioBuffer<float> buffer(2, blockSize);
      juce::MidiBuffer midi;

      setParam(proc, "bypass", 1.0f);
      pumpMessageLoop(20);

      for (int i = 0; i < 50; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }

      auto start = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < numBlocks; ++i) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }
      auto end = std::chrono::high_resolution_clock::now();
      bypassDurationUs =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();

      proc.releaseResources();
    }

    expect(bypassDurationUs > 0, "Bypass duration must be greater than zero");
    expect(activeDurationUs > 0, "Active duration must be greater than zero");

    const double speedupRatio =
        static_cast<double>(activeDurationUs) /
        static_cast<double>(std::max<int64_t>(1, bypassDurationUs));

    expect(speedupRatio >= 2.5,
           "Bypassed processing must be at least 2.5x faster than active "
           "denoising");
  }
};

static PerformanceRatioTest performanceRatioTest;

int main(int argc, char* argv[]) {
  juce::ScopedJuceInitialiser_GUI guiInit;

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
