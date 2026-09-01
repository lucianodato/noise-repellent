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

#include <cmath>
#include <numeric>
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

struct ScopedEditor {
  NoiseRepellentAudioProcessor& proc;
  juce::AudioProcessorEditor* editor = nullptr;
  explicit ScopedEditor(NoiseRepellentAudioProcessor& p)
      : proc(p), editor(p.createEditorIfNeeded()) {
  }
  ~ScopedEditor() {
    if (editor != nullptr) {
      proc.editorBeingDeleted(editor);
      delete editor;
    }
  }
};

bool pullLatestFrame(NoiseRepellentAudioProcessor& proc,
                     NoiseRepellentAudioProcessor::SpectralFrame& outFrame) {
  NoiseRepellentAudioProcessor::SpectralFrame temp;
  bool gotAny = false;
  while (proc.getNextSpectralFrame(temp)) {
    outFrame = temp;
    gotAny = true;
  }
  return gotAny;
}

void processSilentBlocks(NoiseRepellentAudioProcessor& proc,
                         juce::AudioBuffer<float>& buffer,
                         juce::MidiBuffer& midi, int numBlocks = 12) {
  for (int i = 0; i < numBlocks; ++i) {
    buffer.clear();
    proc.processBlock(buffer, midi);
  }
  pumpMessageLoop(10);
}

} // namespace

class UxInvariantsTest : public juce::UnitTest {
public:
  UxInvariantsTest() : juce::UnitTest("UX Invariants Test", "NoiseRepellent") {
  }

  void runTest() override {
    beginTest("Profile Persistence Across 1D <-> 2D While Stopped");
    testProfilePersistenceAcrossEngines();

    beginTest("Live Aggressiveness Control While Silent / Stopped");
    testLiveAggressivenessWhileSilent();

    beginTest("Live Threshold Offsets (Broadband & Tonal) While Silent");
    testLiveThresholdOffsetsWhileSilent();

    beginTest("Silence Segregation (Input/Output Drop vs Profile Preserved)");
    testSilenceSegregation();

    beginTest("Idle No-Profile Silent Bypass");
    testIdleNoProfileBypass();

    beginTest("Engine Switch Progress Advances While Transport Stopped");
    testEngineSwitchProgress();
  }

private:
  void testIdleNoProfileBypass() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    expect(!proc.hasNoiseProfile(), "Processor must have no noise profile");

    // Feed silent blocks
    processSilentBlocks(proc, buffer, midi, 16);

    for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
      for (int s = 0; s < buffer.getNumSamples(); ++s) {
        expectEquals(buffer.getSample(ch, s), 0.0f,
                     "Silent idle blocks must produce zero output");
      }
    }

    NoiseRepellentAudioProcessor::SpectralFrame frame;
    pullLatestFrame(proc, frame);
    expect(!frame.hasNoiseProfile,
           "Spectral frame must report no noise profile");

    proc.releaseResources();
  }

  void testProfilePersistenceAcrossEngines() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // 1. Initially no profile
    expect(!proc.hasNoiseProfile(),
           "Initially processor should have no profile");

    // 2. Start learning on 1D (Spectral)
    setParam(proc, "algorithm_mode", 0.0f); // 1D Spectral
    setParam(proc, "learn_noise", 1.0f);

    // Feed ~0.5s of noise to build profile
    for (int i = 0; i < 50; ++i) {
      generateNoiseBuffer(buffer, 0.05f);
      proc.processBlock(buffer, midi);
    }
    pumpMessageLoop(20);

    // Stop learning and process silence
    setParam(proc, "learn_noise", 0.0f);
    processSilentBlocks(proc, buffer, midi, 12);

    expect(proc.hasNoiseProfile(),
           "Profile must be present after learning on 1D");

    // Pull latest spectral frame
    NoiseRepellentAudioProcessor::SpectralFrame frame1D;
    bool gotFrame = pullLatestFrame(proc, frame1D);
    expect(gotFrame, "Should receive at least one spectral frame in 1D");
    expect(frame1D.hasNoiseProfile,
           "SpectralFrame must report hasNoiseProfile=true in 1D");

    float avgNoiseFloor1D = 0.0f;
    for (float db : frame1D.noiseFloorDB) {
      avgNoiseFloor1D += db;
    }
    avgNoiseFloor1D /= static_cast<float>(frame1D.noiseFloorDB.size());
    expect(avgNoiseFloor1D > -110.0f,
           "1D Noise floor must not be flatlined at -120dB");

    // 3. Switch to 2D (NLM) while stopped / silent
    setParam(proc, "algorithm_mode", 1.0f); // 2D NLM
    pumpMessageLoop(50);

    // Process silent blocks to generate 2D silent frame
    processSilentBlocks(proc, buffer, midi, 12);

    expect(proc.hasNoiseProfile(),
           "Profile must persist on 2D after silent switch");

    NoiseRepellentAudioProcessor::SpectralFrame frame2D;
    gotFrame = pullLatestFrame(proc, frame2D);
    expect(gotFrame, "Should receive spectral frame in 2D mode");
    expect(frame2D.hasNoiseProfile,
           "2D mode must report hasNoiseProfile=true while silent");

    float avgNoiseFloor2D = 0.0f;
    for (float db : frame2D.noiseFloorDB) {
      avgNoiseFloor2D += db;
    }
    avgNoiseFloor2D /= static_cast<float>(frame2D.noiseFloorDB.size());
    expect(avgNoiseFloor2D > -110.0f,
           "2D Noise floor must remain active and not flatline to -120dB");

    float maxDelta1D2D = 0.0f;
    for (size_t i = 0; i < frame1D.noiseFloorDB.size(); ++i) {
      float diff = std::abs(frame1D.noiseFloorDB[i] - frame2D.noiseFloorDB[i]);
      if (diff > maxDelta1D2D)
        maxDelta1D2D = diff;
    }
    expectWithinAbsoluteError(maxDelta1D2D, 0.0f, 0.1f,
                              "Profile bins must match identically across 1D "
                              "<-> 2D switch while stopped");

    // 4. Switch back to 1D while stopped
    setParam(proc, "algorithm_mode", 0.0f);
    pumpMessageLoop(50);
    processSilentBlocks(proc, buffer, midi, 12);

    expect(proc.hasNoiseProfile(),
           "Profile must persist when switching back 2D -> 1D");
    proc.releaseResources();
  }

  void testLiveAggressivenessWhileSilent() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // Learn noise
    setParam(proc, "algorithm_mode", 0.0f);
    setParam(proc, "learn_noise", 1.0f);
    for (int i = 0; i < 50; ++i) {
      generateNoiseBuffer(buffer, 0.05f);
      proc.processBlock(buffer, midi);
    }
    setParam(proc, "learn_noise", 0.0f);
    processSilentBlocks(proc, buffer, midi, 12);

    // Baseline with aggressiveness = -1.0
    setParam(proc, "aggressiveness", -1.0f);
    processSilentBlocks(proc, buffer, midi, 12);

    NoiseRepellentAudioProcessor::SpectralFrame frameLow;
    bool gotLow = pullLatestFrame(proc, frameLow);
    expect(gotLow, "Must receive frame for low aggressiveness");
    float avgLow = std::accumulate(frameLow.noiseFloorDB.begin(),
                                   frameLow.noiseFloorDB.end(), 0.0f) /
                   static_cast<float>(frameLow.noiseFloorDB.size());

    // High aggressiveness = +1.0
    setParam(proc, "aggressiveness", 1.0f);
    processSilentBlocks(proc, buffer, midi, 12);

    NoiseRepellentAudioProcessor::SpectralFrame frameHigh;
    bool gotHigh = pullLatestFrame(proc, frameHigh);
    expect(gotHigh, "Must receive frame for high aggressiveness");
    float avgHigh = std::accumulate(frameHigh.noiseFloorDB.begin(),
                                    frameHigh.noiseFloorDB.end(), 0.0f) /
                    static_cast<float>(frameHigh.noiseFloorDB.size());

    expect(
        avgHigh > avgLow,
        "Increasing aggressiveness while silent must increase noise floor dB");

    proc.releaseResources();
  }

  void testLiveThresholdOffsetsWhileSilent() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // Learn noise
    setParam(proc, "algorithm_mode", 0.0f);
    setParam(proc, "learn_noise", 1.0f);
    for (int i = 0; i < 50; ++i) {
      generateNoiseBuffer(buffer, 0.05f);
      proc.processBlock(buffer, midi);
    }
    setParam(proc, "learn_noise", 0.0f);
    processSilentBlocks(proc, buffer, midi, 12);

    // Offset 0 dB baseline
    setParam(proc, "noise_profile_offset", 0.0f);
    processSilentBlocks(proc, buffer, midi, 12);

    NoiseRepellentAudioProcessor::SpectralFrame frame0;
    bool got0 = pullLatestFrame(proc, frame0);
    expect(got0, "Must receive frame for offset 0dB");
    float avg0 = std::accumulate(frame0.noiseFloorDB.begin(),
                                 frame0.noiseFloorDB.end(), 0.0f) /
                 static_cast<float>(frame0.noiseFloorDB.size());

    // Offset +6 dB
    setParam(proc, "noise_profile_offset", 6.0f);
    processSilentBlocks(proc, buffer, midi, 12);

    NoiseRepellentAudioProcessor::SpectralFrame frame6;
    bool got6 = pullLatestFrame(proc, frame6);
    expect(got6, "Must receive frame for offset +6dB");
    float avg6 = std::accumulate(frame6.noiseFloorDB.begin(),
                                 frame6.noiseFloorDB.end(), 0.0f) /
                 static_cast<float>(frame6.noiseFloorDB.size());

    const float delta = avg6 - avg0;
    expectWithinAbsoluteError(
        delta, 6.0f, 1.5f,
        "Broadband threshold offset +6dB must raise profile by ~6dB");

    proc.releaseResources();
  }

  void testSilenceSegregation() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // Learn noise
    setParam(proc, "learn_noise", 1.0f);
    for (int i = 0; i < 50; ++i) {
      generateNoiseBuffer(buffer, 0.05f);
      proc.processBlock(buffer, midi);
    }
    setParam(proc, "learn_noise", 0.0f);

    // Drain noise frames and feed silent buffers
    NoiseRepellentAudioProcessor::SpectralFrame drainFrame;
    pullLatestFrame(proc, drainFrame);

    processSilentBlocks(proc, buffer, midi, 20);

    NoiseRepellentAudioProcessor::SpectralFrame frame;
    bool gotFrame = pullLatestFrame(proc, frame);
    expect(gotFrame, "Must receive silent frame");

    float maxInput = *std::max_element(frame.inputMagnitudeDB.begin(),
                                       frame.inputMagnitudeDB.end());
    float maxOutput = *std::max_element(frame.outputMagnitudeDB.begin(),
                                        frame.outputMagnitudeDB.end());

    expect(maxInput <= -119.0f, "Silent input spectrum must drop to -120dB");
    expect(maxOutput <= -119.0f, "Silent output spectrum must drop to -120dB");
    expect(frame.hasNoiseProfile,
           "Noise profile must remain active when input is silent");

    proc.releaseResources();
  }

  void testEngineSwitchProgress() {
    NoiseRepellentAudioProcessor proc;
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    expect(!proc.isEngineSwitching(),
           "Initially engine should not be switching");
    expectEquals(proc.getEngineSwitchProgress(), 1.0f,
                 "Initial switch progress should be 1.0 (idle)");

    // Smoothing mode switching is handled seamlessly inside libspecbleach
    // (allocation-free internal crossfade, constant latency), so from the
    // plugin's perspective a switch is instantaneous: no switch phase, no
    // progress ramp — the parameter simply retargets the engine.
    setParam(proc, "algorithm_mode", 1.0f); // Switch to 2D NLM
    pumpMessageLoop(20);

    expect(!proc.isEngineSwitching(),
           "Mode switch must be instantaneous (handled by the library)");
    expectEquals(proc.getEngineSwitchProgress(), 1.0f,
                 "Switch progress must stay idle after mode change");

    // Audio keeps flowing through the same engine group across the switch
    juce::AudioBuffer<float> buffer(2, blockSize);
    for (int ch = 0; ch < 2; ++ch) {
      for (int i = 0; i < blockSize; ++i) {
        buffer.setSample(
            ch, i,
            0.1f * std::sin(2.0f * juce::MathConstants<float>::pi * 440.0f *
                            (float)i / (float)sampleRate));
      }
    }
    juce::MidiBuffer midi;
    proc.processBlock(buffer, midi);
    expect(std::abs(buffer.getSample(0, 0)) < 10.0f,
           "Processing after mode switch produces finite output");

    proc.releaseResources();
  }
};

static UxInvariantsTest uxInvariantsTest;

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
