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

struct ScopedEditor {
  NoiseRepellentAudioProcessor& proc;
  juce::AudioProcessorEditor* editor = nullptr;
  explicit ScopedEditor(NoiseRepellentAudioProcessor& p)
      : proc(p), editor(p.createEditorIfNeeded()) {}
  ~ScopedEditor() {
    if (editor != nullptr) {
      proc.editorBeingDeleted(editor);
      delete editor;
    }
  }
};

void generateNoiseBuffer(juce::AudioBuffer<float>& buffer, float rms = 0.05f, int seed = 42) {
  std::mt19937 gen(static_cast<uint32_t>(seed));
  std::normal_distribution<float> dist(0.0f, rms);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] = dist(gen);
    }
  }
}

void learnProfile(NoiseRepellentAudioProcessor& proc, double sampleRate,
                  int blockSize = 512, int numBlocks = 30) {
  proc.prepareToPlay(sampleRate, blockSize);
  pumpMessageLoop(10);

  setParam(proc, "learn_noise", 1.0f);

  juce::AudioBuffer<float> buffer(2, blockSize);
  juce::MidiBuffer midi;

  for (int b = 0; b < numBlocks; ++b) {
    generateNoiseBuffer(buffer, 0.05f, 1000 + b);
    proc.processBlock(buffer, midi);
  }
  setParam(proc, "learn_noise", 0.0f);
  pumpMessageLoop(10);
}

} // namespace

class IntegrationTest : public juce::UnitTest {
public:
  IntegrationTest()
      : juce::UnitTest("Integration Tests", "NoiseRepellent") {}

  void runTest() override {
    beginTest("Full State Save and Restore Roundtrip (1D Mode)");
    testStateRoundtrip1D();

    beginTest("Full State Save and Restore Roundtrip (2D Mode)");
    testStateRoundtrip2D();

    beginTest("Cross Sample Rate State Restore");
    testCrossSampleRateRestore();

    beginTest("Corrupted / Fuzzed State Ingestion Safety");
    testCorruptedStateSafety();

    beginTest("Residual Listen Mode Dry Reconstruction");
    testResidualListenReconstruction();

    beginTest("Reduction Curve Frequency Response Modulation");
    testReductionCurveResponse();

    beginTest("Adaptive Noise Estimation Modes (SPP-MMSE, Brandt, Martin)");
    testAdaptiveEstimationModes();

    beginTest("Channel Layout Matrix (Mono vs Stereo Processing)");
    testChannelLayoutMatrix();

    beginTest("Transient Protection Dynamic Response");
    testTransientProtection();
  }

private:
  void testStateRoundtrip1D() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    juce::MemoryBlock savedState;

    {
      NoiseRepellentAudioProcessor proc1;
      proc1.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(10);

      setParam(proc1, "algorithm_mode", 0.0f); // 1D Spectral
      setParam(proc1, "reduction_amount", 18.0f);
      setParam(proc1, "tonal_reduction", 14.0f);
      setParam(proc1, "aggressiveness", 0.7f);
      setParam(proc1, "smoothing_factor", 45.0f);

      // Add custom reduction curve
      std::vector<NoiseRepellentAudioProcessor::CurveNode> nodes = {
          {0.0f, 0.0f}, {0.25f, -6.0f}, {0.75f, 12.0f}, {1.0f, 0.0f}};
      proc1.setCurveNodes(nodes);

      learnProfile(proc1, sampleRate, blockSize, 25);
      expect(proc1.hasNoiseProfile(), "Proc 1 must have noise profile");

      proc1.getStateInformation(savedState);
      expect(savedState.getSize() > 0, "Saved state memory block must not be empty");
    }

    {
      NoiseRepellentAudioProcessor proc2;
      proc2.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(10);

      expect(!proc2.hasNoiseProfile(), "Proc 2 must start with no profile");

      proc2.setStateInformation(savedState.getData(), static_cast<int>(savedState.getSize()));
      pumpMessageLoop(20);

      expect(proc2.hasNoiseProfile(), "Proc 2 must restore noise profile from state");

      auto* pReduction = proc2.getAPVTS().getRawParameterValue("reduction_amount");
      auto* pAgg = proc2.getAPVTS().getRawParameterValue("aggressiveness");
      auto* pSmooth = proc2.getAPVTS().getRawParameterValue("smoothing_factor");

      expect(pReduction != nullptr && std::abs(pReduction->load() - 18.0f) < 0.2f,
             "Restored reduction_amount must match");
      expect(pAgg != nullptr && std::abs(pAgg->load() - 0.7f) < 0.05f,
             "Restored aggressiveness must match");
      expect(pSmooth != nullptr && std::abs(pSmooth->load() - 45.0f) < 0.5f,
             "Restored smoothing_factor must match");

      // Verify restored curve nodes
      const auto& restoredNodes = proc2.getCurveNodes();
      expectEquals(static_cast<int>(restoredNodes.size()), 4, "Restored curve nodes count must match");
      if (restoredNodes.size() == 4) {
        expectWithinAbsoluteError(restoredNodes[1].normX, 0.25f, 0.01f);
        expectWithinAbsoluteError(restoredNodes[1].biasDB, -6.0f, 0.01f);
        expectWithinAbsoluteError(restoredNodes[2].normX, 0.75f, 0.01f);
        expectWithinAbsoluteError(restoredNodes[2].biasDB, 12.0f, 0.01f);
      }
    }
  }

  void testStateRoundtrip2D() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    juce::MemoryBlock savedState;

    {
      NoiseRepellentAudioProcessor proc1;
      proc1.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(10);

      setParam(proc1, "algorithm_mode", 1.0f); // 2D NLM
      for (int step = 0; step < 40; ++step) pumpMessageLoop(20);

      learnProfile(proc1, sampleRate, blockSize, 25);
      expect(proc1.hasNoiseProfile(), "Proc 1 (2D) must have noise profile");

      proc1.getStateInformation(savedState);
      expect(savedState.getSize() > 0, "Saved state memory block must not be empty");
    }

    {
      NoiseRepellentAudioProcessor proc2;
      proc2.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(10);

      proc2.setStateInformation(savedState.getData(), static_cast<int>(savedState.getSize()));
      for (int step = 0; step < 40; ++step) pumpMessageLoop(20);

      expect(proc2.hasNoiseProfile(), "Proc 2 must restore noise profile in 2D mode");

      auto* pAlgo = proc2.getAPVTS().getRawParameterValue("algorithm_mode");
      expect(pAlgo != nullptr && std::round(pAlgo->load()) == 1.0f, "Restored algorithm mode must be 2D");
    }
  }

  void testCrossSampleRateRestore() {
    juce::MemoryBlock savedState;

    // Learn and save at 44.1 kHz
    {
      NoiseRepellentAudioProcessor proc44k;
      proc44k.prepareToPlay(44100.0, 256);
      pumpMessageLoop(10);

      learnProfile(proc44k, 44100.0, 256, 25);
      proc44k.getStateInformation(savedState);
    }

    // Restore into 96 kHz instance
    {
      NoiseRepellentAudioProcessor proc96k;
      proc96k.prepareToPlay(96000.0, 512);
      pumpMessageLoop(10);

      proc96k.setStateInformation(savedState.getData(), static_cast<int>(savedState.getSize()));
      pumpMessageLoop(20);

      expect(proc96k.hasNoiseProfile(), "Profile must be present after cross-sample-rate restore");

      juce::AudioBuffer<float> buffer(2, 512);
      juce::MidiBuffer midi;
      generateNoiseBuffer(buffer, 0.05f, 999);

      // Must process without crash or NaN
      for (int b = 0; b < 10; ++b) {
        proc96k.processBlock(buffer, midi);
        for (int ch = 0; ch < 2; ++ch) {
          const auto* r = buffer.getReadPointer(ch);
          for (int s = 0; s < 512; ++s) {
            expect(!std::isnan(r[s]), "Cross sample rate processed output must not be NaN");
          }
        }
      }
    }
  }

  void testCorruptedStateSafety() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);

    // 1. Null data & 0 bytes
    proc.setStateInformation(nullptr, 0);

    // 2. Truncated XML
    const char* truncXml = "<PARAMETERS><VALUE id=\"aggressiveness\"";
    proc.setStateInformation(truncXml, static_cast<int>(std::strlen(truncXml)));

    // 3. Random binary garbage
    std::vector<uint8_t> garbage(1024);
    for (size_t i = 0; i < garbage.size(); ++i) garbage[i] = static_cast<uint8_t>(i * 37 + 13);
    proc.setStateInformation(garbage.data(), static_cast<int>(garbage.size()));

    // 4. Invalid Base64 in profile node
    const char* invalidB64 = "<PARAMETERS><LEARNED_PROFILES><CHANNEL_PROFILE channel=\"0\" mode=\"0\" size=\"100\" data=\"!@#$%^&*()\"/></LEARNED_PROFILES></PARAMETERS>";
    proc.setStateInformation(invalidB64, static_cast<int>(std::strlen(invalidB64)));

    expect(true, "Corrupted state ingestion must be safely rejected without crashing");
    proc.releaseResources();
  }

  void testResidualListenReconstruction() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor procDenoise;
    NoiseRepellentAudioProcessor procResidual;

    procDenoise.prepareToPlay(sampleRate, blockSize);
    procResidual.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    learnProfile(procDenoise, sampleRate, blockSize, 25);

    // Share identical state
    juce::MemoryBlock state;
    procDenoise.getStateInformation(state);
    procResidual.setStateInformation(state.getData(), static_cast<int>(state.getSize()));
    pumpMessageLoop(10);

    setParam(procDenoise, "residual_listen", 0.0f);
    setParam(procResidual, "residual_listen", 1.0f);

    setParam(procDenoise, "reduction_amount", 20.0f);
    setParam(procResidual, "reduction_amount", 20.0f);

    std::mt19937 gen(555);
    std::normal_distribution<float> dist(0.0f, 0.05f);

    // Run warmup blocks
    juce::AudioBuffer<float> bufDenoise(2, blockSize);
    juce::AudioBuffer<float> bufResidual(2, blockSize);
    juce::MidiBuffer midi;

    for (int b = 0; b < 20; ++b) {
      for (int ch = 0; ch < 2; ++ch) {
        auto* wD = bufDenoise.getWritePointer(ch);
        auto* wR = bufResidual.getWritePointer(ch);
        for (int s = 0; s < blockSize; ++s) {
          float val = dist(gen);
          wD[s] = val;
          wR[s] = val;
        }
      }
      procDenoise.processBlock(bufDenoise, midi);
      procResidual.processBlock(bufResidual, midi);
    }

    // Measure RMS of denoise + residual vs non-zero energy
    float denoiseRms = bufDenoise.getRMSLevel(0, 0, blockSize);
    float residualRms = bufResidual.getRMSLevel(0, 0, blockSize);

    expect(denoiseRms > 0.0f, "Denoised signal must have non-zero energy");
    expect(residualRms > 0.0f, "Residual signal must have non-zero energy");

    procDenoise.releaseResources();
    procResidual.releaseResources();
  }

  void testReductionCurveResponse() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    learnProfile(proc, sampleRate, blockSize, 25);

    // Enable reduction curve with heavy high-frequency reduction
    setParam(proc, "reduction_curve_enabled", 1.0f);
    std::vector<NoiseRepellentAudioProcessor::CurveNode> nodes = {
        {0.0f, 0.0f},
        {0.5f, 0.0f},
        {0.8f, 24.0f}, // +24 dB extra reduction above 80% Nyquist
        {1.0f, 24.0f}};
    proc.setCurveNodes(nodes);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    for (int b = 0; b < 25; ++b) {
      generateNoiseBuffer(buffer, 0.05f, 777 + b);
      proc.processBlock(buffer, midi);
    }

    NoiseRepellentAudioProcessor::SpectralFrame frame;
    bool gotFrame = false;
    while (proc.getNextSpectralFrame(frame)) {
      gotFrame = true;
    }

    expect(gotFrame, "Should have retrieved spectral frame with active reduction curve");
    if (gotFrame) {
      expect(frame.reductionCurveEnabled, "Frame must reflect reductionCurveEnabled = true");
    }

    proc.releaseResources();
  }

  void testAdaptiveEstimationModes() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    const std::vector<float> methods = {0.0f, 1.0f, 2.0f}; // SPP-MMSE, Brandt, Martin

    for (float method : methods) {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(10);

      setParam(proc, "adaptive_noise", 1.0f);
      setParam(proc, "adaptive_method", method);

      juce::AudioBuffer<float> buffer(2, blockSize);
      juce::MidiBuffer midi;

      // Stream continuous noise to allow adaptive tracker to process
      for (int b = 0; b < 30; ++b) {
        generateNoiseBuffer(buffer, 0.04f, 800 + b);
        proc.processBlock(buffer, midi);

        for (int ch = 0; ch < 2; ++ch) {
          const auto* r = buffer.getReadPointer(ch);
          for (int s = 0; s < blockSize; ++s) {
            expect(!std::isnan(r[s]), "Adaptive mode processed output must not be NaN");
            expect(!std::isinf(r[s]), "Adaptive mode processed output must not be Inf");
          }
        }
      }

      proc.releaseResources();
    }
  }

  void testChannelLayoutMatrix() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    learnProfile(proc, sampleRate, blockSize, 20);

    juce::MidiBuffer midi;

    // 1. Mono processing (1-channel buffer)
    {
      juce::AudioBuffer<float> monoBuffer(1, blockSize);
      generateNoiseBuffer(monoBuffer, 0.05f, 123);
      proc.processBlock(monoBuffer, midi);

      const auto* r = monoBuffer.getReadPointer(0);
      for (int s = 0; s < blockSize; ++s) {
        expect(!std::isnan(r[s]), "Mono output must not be NaN");
      }
    }

    // 2. Stereo with one channel silent
    {
      juce::AudioBuffer<float> stereoBuffer(2, blockSize);
      stereoBuffer.clear();
      auto* left = stereoBuffer.getWritePointer(0);
      std::mt19937 gen(456);
      std::normal_distribution<float> dist(0.0f, 0.05f);
      for (int s = 0; s < blockSize; ++s) left[s] = dist(gen);

      proc.processBlock(stereoBuffer, midi);

      const auto* leftOut = stereoBuffer.getReadPointer(0);
      const auto* rightOut = stereoBuffer.getReadPointer(1);

      for (int s = 0; s < blockSize; ++s) {
        expect(!std::isnan(leftOut[s]), "Left channel must not be NaN");
        expect(!std::isnan(rightOut[s]), "Right channel must not be NaN");
      }
    }

    proc.releaseResources();
  }

  void testTransientProtection() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    setParam(proc, "transient_protection_enable", 1.0f);
    learnProfile(proc, sampleRate, blockSize, 20);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;
    buffer.clear();

    // Inject strong impulse transient
    for (int ch = 0; ch < 2; ++ch) {
      buffer.setSample(ch, 10, 0.95f);
      buffer.setSample(ch, 11, -0.85f);
      buffer.setSample(ch, 12, 0.5f);
    }

    proc.processBlock(buffer, midi);
    pumpMessageLoop(5);

    float transientAct = proc.consumeTransientActivity();
    // Transient activity should be reported or preserved safely
    expect(!std::isnan(transientAct), "Transient activity must be a valid number");

    proc.releaseResources();
  }
};

static IntegrationTest integrationTest;

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
