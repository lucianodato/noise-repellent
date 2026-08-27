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

#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <thread>
#include <vector>

#include <juce_core/juce_core.h>
#include <juce_events/juce_events.h>
#include <juce_gui_basics/juce_gui_basics.h>

#include "PluginProcessor.h"

namespace {

void pumpMessageLoop(int maxMillis = 20) {
  if (auto* mm = juce::MessageManager::getInstanceWithoutCreating()) {
    mm->runDispatchLoopUntil(maxMillis);
  }
}

void setParam(NoiseRepellentAudioProcessor& proc, const juce::String& paramId,
              float value) {
  if (auto* param = proc.getAPVTS().getParameter(paramId)) {
    param->setValueNotifyingHost(param->convertTo0to1(value));
  }
  pumpMessageLoop(2);
}

void learnProfile(NoiseRepellentAudioProcessor& proc, double sampleRate,
                  int blockSize = 512, int numBlocks = 30) {
  proc.prepareToPlay(sampleRate, blockSize);
  pumpMessageLoop(10);

  setParam(proc, "learn_noise", 1.0f);

  std::mt19937 gen(42);
  std::normal_distribution<float> dist(0.0f, 0.05f);

  juce::AudioBuffer<float> buffer(2, blockSize);
  juce::MidiBuffer midi;

  for (int b = 0; b < numBlocks; ++b) {
    for (int ch = 0; ch < 2; ++ch) {
      auto* data = buffer.getWritePointer(ch);
      for (int s = 0; s < blockSize; ++s) {
        data[s] = dist(gen);
      }
    }
    proc.processBlock(buffer, midi);
  }
  setParam(proc, "learn_noise", 0.0f);
  pumpMessageLoop(10);
}

} // namespace

class StressAndRobustnessTest : public juce::UnitTest {
public:
  StressAndRobustnessTest()
      : juce::UnitTest("Stress & Robustness Tests", "NoiseRepellent") {}

  void runTest() override {
    beginTest("Sample Rate Switching Fuzzing Under Continuous Audio Load");
    testSampleRateSwitchingStress();

    beginTest("Dynamic & Arbitrary Non-Power-of-2 Buffer Size Fuzzing");
    testArbitraryBlockSizeStress();

    beginTest("Pathological Signal Robustness (Denormals, DC Offset, Clipping, Bursts, Silence)");
    testPathologicalSignalTorture();

    beginTest("Multithreaded Concurrency & Parameter Storm");
    testMultithreadedParameterStorm();

    beginTest("Rapid Prepare/Release/Reset Lifecycle Thrashing");
    testPrepareReleaseThrashing();
  }

private:
  void testSampleRateSwitchingStress() {
    NoiseRepellentAudioProcessor proc;
    const std::vector<double> sampleRates = {44100.0, 48000.0, 88200.0, 96000.0, 192000.0, 44100.0};
    const std::vector<int> blockSizes = {64, 128, 256, 512, 1024};

    std::mt19937 gen(101);
    std::normal_distribution<float> dist(0.0f, 0.03f);
    juce::MidiBuffer midi;

    learnProfile(proc, 48000.0, 512, 20);
    expect(proc.hasNoiseProfile(), "Profile must be present after learning");

    for (int iter = 0; iter < 12; ++iter) {
      double sr = sampleRates[static_cast<size_t>(iter % sampleRates.size())];
      int bs = blockSizes[static_cast<size_t>(iter % blockSizes.size())];

      proc.prepareToPlay(sr, bs);
      pumpMessageLoop(5);

      juce::AudioBuffer<float> buffer(2, bs);
      for (int b = 0; b < 10; ++b) {
        for (int ch = 0; ch < 2; ++ch) {
          auto* w = buffer.getWritePointer(ch);
          for (int s = 0; s < bs; ++s) {
            w[s] = dist(gen);
          }
        }
        proc.processBlock(buffer, midi);

        for (int ch = 0; ch < 2; ++ch) {
          const auto* r = buffer.getReadPointer(ch);
          for (int s = 0; s < bs; ++s) {
            expect(!std::isnan(r[s]), "Output sample must not be NaN during sample rate switch");
            expect(!std::isinf(r[s]), "Output sample must not be Inf during sample rate switch");
          }
        }
      }
    }
    proc.releaseResources();
  }

  void testArbitraryBlockSizeStress() {
    NoiseRepellentAudioProcessor proc;
    const double sr = 48000.0;
    proc.prepareToPlay(sr, 4096);
    pumpMessageLoop(10);

    learnProfile(proc, sr, 512, 20);

    // Arbitrary, prime, non-power-of-2 block sizes
    const std::vector<int> oddBlockSizes = {
        1, 7, 13, 17, 31, 63, 100, 255, 333, 513, 777, 1023, 1500, 2049, 3999, 4096
    };

    std::mt19937 gen(202);
    std::normal_distribution<float> dist(0.0f, 0.04f);
    juce::MidiBuffer midi;

    for (int bs : oddBlockSizes) {
      juce::AudioBuffer<float> buffer(2, bs);
      for (int b = 0; b < 6; ++b) {
        for (int ch = 0; ch < 2; ++ch) {
          auto* w = buffer.getWritePointer(ch);
          for (int s = 0; s < bs; ++s) {
            w[s] = dist(gen);
          }
        }
        proc.processBlock(buffer, midi);

        for (int ch = 0; ch < 2; ++ch) {
          const auto* r = buffer.getReadPointer(ch);
          for (int s = 0; s < bs; ++s) {
            expect(!std::isnan(r[s]), "Output sample must not be NaN for block size " + juce::String(bs));
            expect(!std::isinf(r[s]), "Output sample must not be Inf for block size " + juce::String(bs));
          }
        }
      }
    }
    proc.releaseResources();
  }

  void testPathologicalSignalTorture() {
    constexpr double sr = 48000.0;
    constexpr int bs = 256;

    auto testSignal = [&](const juce::String& description, auto generator) {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(sr, bs);
      pumpMessageLoop(10);
      learnProfile(proc, sr, bs, 15);

      juce::AudioBuffer<float> buffer(2, bs);
      juce::MidiBuffer midi;

      for (int b = 0; b < 10; ++b) {
        for (int ch = 0; ch < 2; ++ch) {
          auto* w = buffer.getWritePointer(ch);
          for (int s = 0; s < bs; ++s) {
            w[s] = generator(ch, s, b);
          }
        }
        proc.processBlock(buffer, midi);

        for (int ch = 0; ch < 2; ++ch) {
          const auto* r = buffer.getReadPointer(ch);
          for (int s = 0; s < bs; ++s) {
            expect(!std::isnan(r[s]), description + ": output must not be NaN");
            expect(!std::isinf(r[s]), description + ": output must not be Inf");
          }
        }
      }
      proc.releaseResources();
    };

    // 1. Subnormals / Denormals
    testSignal("Denormal input", [](int, int, int) {
      return std::numeric_limits<float>::denorm_min();
    });

    // 2. Extreme DC Offset (+100.0f and -100.0f)
    testSignal("Extreme DC offset", [](int ch, int, int) {
      return ch == 0 ? 100.0f : -100.0f;
    });

    // 3. Full scale square wave (+50 dBFS clipping burst)
    testSignal("Extreme full-scale clipping", [](int, int s, int) {
      return (s % 4 < 2) ? 500.0f : -500.0f;
    });

    // 4. Nyquist alternating sequence (+1.0, -1.0, +1.0, -1.0)
    testSignal("Nyquist tone", [](int, int s, int) {
      return (s % 2 == 0) ? 0.9f : -0.9f;
    });

    // 5. Complete Digital Silence
    testSignal("Digital silence", [](int, int, int) {
      return 0.0f;
    });

    // 6. Crash resistance on NaNs and Infs
    {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(sr, bs);
      juce::AudioBuffer<float> nanBuf(2, bs);
      juce::MidiBuffer midi;
      for (int ch = 0; ch < 2; ++ch) {
        auto* w = nanBuf.getWritePointer(ch);
        for (int s = 0; s < bs; ++s) {
          w[s] = (s % 2 == 0) ? std::numeric_limits<float>::quiet_NaN() : std::numeric_limits<float>::infinity();
        }
      }
      // Must not crash or throw unhandled exception
      proc.processBlock(nanBuf, midi);
      proc.releaseResources();
      expect(true, "NaN / Inf ingestion must not cause processor crash");
    }
  }

  void testMultithreadedParameterStorm() {
    NoiseRepellentAudioProcessor proc;
    const double sr = 48000.0;
    const int bs = 256;
    proc.prepareToPlay(sr, bs);
    pumpMessageLoop(10);

    learnProfile(proc, sr, bs, 20);

    std::atomic<bool> keepRunning{true};
    std::atomic<int> blocksProcessed{0};

    // Audio rendering thread
    std::thread audioThread([&]() {
      juce::AudioBuffer<float> buffer(2, bs);
      juce::MidiBuffer midi;
      std::mt19937 gen(303);
      std::normal_distribution<float> dist(0.0f, 0.05f);

      while (keepRunning.load(std::memory_order_relaxed)) {
        for (int ch = 0; ch < 2; ++ch) {
          auto* w = buffer.getWritePointer(ch);
          for (int s = 0; s < bs; ++s) {
            w[s] = dist(gen);
          }
        }
        proc.processBlock(buffer, midi);
        blocksProcessed.fetch_add(1, std::memory_order_relaxed);
      }
    });

    // Parameter storm on host/GUI thread
    std::mt19937 paramGen(404);
    std::uniform_real_distribution<float> unitDist(0.0f, 1.0f);
    std::uniform_real_distribution<float> aggDist(-1.0f, 1.0f);
    std::uniform_real_distribution<float> redDist(0.0f, 40.0f);
    std::uniform_real_distribution<float> offDist(-12.0f, 12.0f);

    const auto startTime = std::chrono::steady_clock::now();
    while (std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now() - startTime)
               .count() < 1200) {
      setParam(proc, "aggressiveness", aggDist(paramGen));
      setParam(proc, "reduction_amount", redDist(paramGen));
      setParam(proc, "tonal_reduction", redDist(paramGen));
      setParam(proc, "noise_profile_offset", offDist(paramGen));
      setParam(proc, "tonal_noise_profile_offset", offDist(paramGen));
      setParam(proc, "smoothing_factor", unitDist(paramGen) * 100.0f);
      setParam(proc, "masking_depth", unitDist(paramGen) * 100.0f);
      setParam(proc, "whitening_factor", unitDist(paramGen) * 100.0f);
      setParam(proc, "transient_protection_enable", unitDist(paramGen) > 0.5f ? 1.0f : 0.0f);
      setParam(proc, "reduction_curve_enabled", unitDist(paramGen) > 0.5f ? 1.0f : 0.0f);

      // Random curve node updates
      if (unitDist(paramGen) > 0.7f) {
        std::vector<NoiseRepellentAudioProcessor::CurveNode> nodes = {
            {0.0f, offDist(paramGen)},
            {0.3f, offDist(paramGen)},
            {0.7f, offDist(paramGen)},
            {1.0f, offDist(paramGen)}};
        proc.setCurveNodes(nodes);
      }

      pumpMessageLoop(1);
    }

    keepRunning.store(false, std::memory_order_relaxed);
    audioThread.join();

    expect(blocksProcessed.load() > 100, "Audio thread must have processed blocks during parameter storm");
    proc.releaseResources();
  }

  void testPrepareReleaseThrashing() {
    NoiseRepellentAudioProcessor proc;
    const double sr = 48000.0;
    const int bs = 512;

    for (int i = 0; i < 20; ++i) {
      proc.prepareToPlay(sr, bs);
      pumpMessageLoop(2);

      juce::AudioBuffer<float> buffer(2, bs);
      juce::MidiBuffer midi;
      buffer.clear();
      proc.processBlock(buffer, midi);

      proc.reset();
      proc.releaseResources();
      pumpMessageLoop(2);
    }
    expect(true, "Prepare/Release/Reset thrashing must complete without crash");
  }
};

static StressAndRobustnessTest stressAndRobustnessTest;

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
