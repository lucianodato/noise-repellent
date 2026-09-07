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

void learnNoiseProfile(NoiseRepellentAudioProcessor& proc,
                       juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midi,
                       int numBlocks = 50) {
  setParam(proc, "learn_noise", 1.0f);
  for (int i = 0; i < numBlocks; ++i) {
    generateNoiseBuffer(buffer, 0.05f);
    proc.processBlock(buffer, midi);
  }
  setParam(proc, "learn_noise", 0.0f);
  pumpMessageLoop(20);
}

void generateTonalBuffer(juce::AudioBuffer<float>& buffer, double sampleRate,
                         float freqHz = 440.0f, float amplitude = 0.1f) {
  static double phase = 0.0;
  const double inc = 2.0 * 3.141592653589793 * freqHz / sampleRate;
  std::mt19937 gen(4242);
  std::normal_distribution<float> dist(0.0f, 0.01f);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] =
          amplitude * static_cast<float>(std::sin(phase)) + dist(gen);
      phase += inc;
    }
  }
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

    beginTest("Low-Latency Enter/Exit Drops Profile And Locks Controls");
    testLowLatencyEnterExitCleanSlate();

    beginTest("Frame-Size Switch Clean Slate And Learn Auto-Stop");
    testFrameSizeSwitchCleanSlate();

    beginTest("Tonal Peak Visibility Rule And State Preservation");
    testTonalPeakVisibilityRule();

    beginTest("Mode Switch Never Suspends Or Re-Reports Latency");
    testModeSwitchNeverSuspends();
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

    // 3b. Switch to NLM+DFTT while stopped: intra-family flip, profile and
    // floor must survive identically
    setParam(proc, "algorithm_mode", 2.0f); // NLM + DFTT refinement
    pumpMessageLoop(50);
    processSilentBlocks(proc, buffer, midi, 12);

    expect(proc.hasNoiseProfile(),
           "Profile must persist on NLM+DFTT after silent switch");

    NoiseRepellentAudioProcessor::SpectralFrame frameDFTT;
    gotFrame = pullLatestFrame(proc, frameDFTT);
    expect(gotFrame, "Should receive spectral frame in NLM+DFTT mode");
    expect(frameDFTT.hasNoiseProfile,
           "NLM+DFTT mode must report hasNoiseProfile=true while silent");

    float avgNoiseFloorDFTT = 0.0f;
    for (float db : frameDFTT.noiseFloorDB) {
      avgNoiseFloorDFTT += db;
    }
    avgNoiseFloorDFTT /= static_cast<float>(frameDFTT.noiseFloorDB.size());
    expect(avgNoiseFloorDFTT > -110.0f,
           "NLM+DFTT noise floor must remain active and not flatline");

    float maxDelta2DDFTT = 0.0f;
    for (size_t i = 0; i < frame2D.noiseFloorDB.size(); ++i) {
      float diff =
          std::abs(frame2D.noiseFloorDB[i] - frameDFTT.noiseFloorDB[i]);
      if (diff > maxDelta2DDFTT)
        maxDelta2DDFTT = diff;
    }
    expectWithinAbsoluteError(maxDelta2DDFTT, 0.0f, 0.1f,
                              "Profile bins must match identically across 2D "
                              "<-> NLM+DFTT switch while stopped");

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

  void testLowLatencyEnterExitCleanSlate() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // Unlink first so the forced-link behavior is observable
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    setParam(proc, "algorithm_mode", 0.0f);
    learnNoiseProfile(proc, buffer, midi);
    processSilentBlocks(proc, buffer, midi, 12);
    expect(proc.hasNoiseProfile(), "Profile must exist before low-latency");
    const int normalLatency = proc.getLatencySamples();
    expect(normalLatency > 512, "Normal latency must exceed one 512 frame");

    // Enter low-latency: clean slate + defaults + re-reported PDC
    setParam(proc, "low_latency", 1.0f);
    pumpMessageLoop(50);
    processSilentBlocks(proc, buffer, midi, 12);

    expect(proc.isLowLatency(), "Processor must report low-latency active");
    expect(!proc.hasNoiseProfile(),
           "Entering low-latency must drop the learned profile");
    expect(proc.getAPVTS().getRawParameterValue("learn_noise")->load() < 0.5f,
           "Entering low-latency must auto-stop Learn");
    expectEquals(proc.getLatencySamples(), 512,
                 "Low-latency must re-report one 512-sample frame");

    NoiseRepellentAudioProcessor::SpectralFrame frame;
    pullLatestFrame(proc, frame);
    expect(frame.isLinked && frame.isOffsetLinked,
           "Low-latency must force reduction/threshold links on");

    // Exit: still clean, latency restored, links follow parameters again
    setParam(proc, "low_latency", 0.0f);
    pumpMessageLoop(50);
    processSilentBlocks(proc, buffer, midi, 12);

    expect(!proc.isLowLatency(), "Low-latency must be off after exit");
    expect(!proc.hasNoiseProfile(),
           "Exiting low-latency must not resurrect a profile");
    expectEquals(proc.getLatencySamples(), normalLatency,
                 "Exiting low-latency must restore normal latency");
    pullLatestFrame(proc, frame);
    expect(frame.isLinked && frame.isOffsetLinked,
           "After exit, forced links persist until the user unlinks again");

    proc.releaseResources();
  }

  void testFrameSizeSwitchCleanSlate() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // Learn, then leave Learn running: the switch must auto-stop it
    setParam(proc, "learn_noise", 1.0f);
    for (int i = 0; i < 50; ++i) {
      generateNoiseBuffer(buffer, 0.05f);
      proc.processBlock(buffer, midi);
    }
    pumpMessageLoop(20);

    auto* choice = dynamic_cast<juce::AudioParameterChoice*>(
        proc.getAPVTS().getParameter("frame_size"));
    expect(choice != nullptr, "frame_size choice parameter must exist");
    const float otherIndex = choice->getIndex() == 0 ? 4.0f : 0.0f;
    const int latencyBefore = proc.getLatencySamples();
    setParam(proc, "frame_size", otherIndex);
    pumpMessageLoop(50);
    processSilentBlocks(proc, buffer, midi, 12);

    expect(!proc.hasNoiseProfile(),
           "Frame-size switch must discard the learned profile");
    expect(proc.getAPVTS().getRawParameterValue("learn_noise")->load() < 0.5f,
           "Frame-size switch must auto-stop an in-progress Learn");
    expect(proc.getLatencySamples() != latencyBefore,
           "Frame-size switch must re-report latency via suspended rebuild");

    NoiseRepellentAudioProcessor::SpectralFrame frame;
    pullLatestFrame(proc, frame);
    expect(!frame.hasNoiseProfile,
           "HUD frame must show NO PROFILE after frame-size switch");

    proc.releaseResources();
  }

  void testTonalPeakVisibilityRule() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // Learn on tonal material so peaks exist; links on by default
    setParam(proc, "link_reduction", 1.0f);
    setParam(proc, "link_threshold_offset", 1.0f);
    setParam(proc, "learn_noise", 1.0f);
    for (int i = 0; i < 50; ++i) {
      generateTonalBuffer(buffer, sampleRate);
      proc.processBlock(buffer, midi);
    }
    setParam(proc, "learn_noise", 0.0f);
    processSilentBlocks(proc, buffer, midi, 12);
    expect(proc.hasNoiseProfile(), "Tonal profile must be learned");

    // Linked + profile: peaks hidden (visibility rule, "only if" direction)
    NoiseRepellentAudioProcessor::SpectralFrame linkedFrame;
    pullLatestFrame(proc, linkedFrame);
    expect(linkedFrame.tonalPeaksHz.empty(),
           "Tonal peaks must be hidden while fully linked");

    // Unlinked + profile: peaks visible
    setParam(proc, "link_threshold_offset", 0.0f);
    processSilentBlocks(proc, buffer, midi, 12);
    NoiseRepellentAudioProcessor::SpectralFrame unlinkedFrame;
    pullLatestFrame(proc, unlinkedFrame);
    expect(!unlinkedFrame.tonalPeaksHz.empty(),
           "Tonal peaks must be visible once unlinked with a profile");

    // State preservation: identical peak set across a mode switch + silence
    setParam(proc, "algorithm_mode", 1.0f);
    pumpMessageLoop(50);
    processSilentBlocks(proc, buffer, midi, 12);
    NoiseRepellentAudioProcessor::SpectralFrame switchedFrame;
    pullLatestFrame(proc, switchedFrame);
    expect(!switchedFrame.tonalPeaksHz.empty(),
           "Tonal peaks must survive smoothing mode switches");
    expectEquals(switchedFrame.tonalPeaksHz.size(),
                 unlinkedFrame.tonalPeaksHz.size(),
                 "Peak set must be preserved across mode switch");

    proc.releaseResources();
  }

  void testModeSwitchNeverSuspends() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    setParam(proc, "algorithm_mode", 0.0f);
    learnNoiseProfile(proc, buffer, midi);
    processSilentBlocks(proc, buffer, midi, 12);
    const int latency = proc.getLatencySamples();

    for (float mode : {1.0f, 2.0f, 0.0f}) {
      setParam(proc, "algorithm_mode", mode);
      pumpMessageLoop(50);
      processSilentBlocks(proc, buffer, midi, 12);

      expect(!proc.isSuspended(),
             "Smoothing mode switch must never suspend processing");
      expectEquals(proc.getLatencySamples(), latency,
                   "Smoothing mode switch must not re-report latency");
      expect(proc.hasNoiseProfile(),
             "Profile must survive smoothing mode switches");
    }

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
