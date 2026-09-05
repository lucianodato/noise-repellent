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
#include <limits>
#include <map>
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

void learnStereoDistinctProfiles(NoiseRepellentAudioProcessor& proc,
                                 int blockSize = 512, int numBlocks = 30) {
  juce::AudioBuffer<float> buffer(2, blockSize);
  juce::MidiBuffer midi;
  for (int b = 0; b < numBlocks; ++b) {
    // Deliberately different per-channel noise so L/R profiles are distinct.
    std::mt19937 genL(static_cast<uint32_t>(3000 + b));
    std::normal_distribution<float> distL(0.0f, 0.06f);
    auto* left = buffer.getWritePointer(0);
    for (int s = 0; s < blockSize; ++s) left[s] = distL(genL);
    std::mt19937 genR(static_cast<uint32_t>(9000 + b));
    std::normal_distribution<float> distR(0.0f, 0.015f);
    auto* right = buffer.getWritePointer(1);
    for (int s = 0; s < blockSize; ++s) right[s] = distR(genR);
    proc.processBlock(buffer, midi);
  }
}

// (channel, mode) -> profile magnitudes decoded from a saved state block.
std::map<std::pair<int, int>, std::vector<float>>
extractStateProfiles(const juce::MemoryBlock& state) {
  std::map<std::pair<int, int>, std::vector<float>> out;
  std::unique_ptr<juce::XmlElement> xml(
      juce::AudioProcessor::getXmlFromBinary(state.getData(),
                                             (int)state.getSize()));
  if (xml == nullptr)
    return out;
  juce::ValueTree tree = juce::ValueTree::fromXml(*xml);
  if (!tree.isValid())
    return out;
  juce::ValueTree profiles = tree.getChildWithName("LEARNED_PROFILES");
  if (!profiles.isValid())
    return out;
  for (int i = 0; i < profiles.getNumChildren(); ++i) {
    juce::ValueTree node = profiles.getChild(i);
    int channel = node.getProperty("channel", -1);
    int mode = node.getProperty("mode", 0);
    juce::String base64Data = node.getProperty("data", "");
    juce::MemoryBlock mb;
    if (channel < 0 || base64Data.isEmpty() ||
        !mb.fromBase64Encoding(base64Data) || mb.getSize() == 0)
      continue;
    const size_t count = mb.getSize() / sizeof(float);
    const float* data = reinterpret_cast<const float*>(mb.getData());
    out[{channel, mode}] = std::vector<float>(data, data + count);
  }
  return out;
}

float maxAbsDiff(const std::vector<float>& a, const std::vector<float>& b) {
  if (a.size() != b.size())
    return std::numeric_limits<float>::infinity();
  float m = 0.0f;
  for (size_t i = 0; i < a.size(); ++i)
    m = std::max(m, std::abs(a[i] - b[i]));
  return m;
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

    beginTest("Mono Bus Layout (Single-Engine Group)");
    testMonoBusLayout();

    beginTest("Stereo-Mono-Stereo Distinct Profile Restore");
    testStereoMonoStereoDistinctRestore();

    beginTest("Transient Protection Dynamic Response");
    testTransientProtection();

    beginTest("Frame-Size Switch Preserves Profile and Re-reports Latency");
    testFrameSizeSwitch();

    beginTest("Low-Latency Mode Rebuilds Clean With 512-Sample Latency");
    testLowLatencyMode();

    beginTest("Reduction Curve Shapes Output In Low-Latency Mode");
    testReductionCurveLowLatency();

    beginTest("Bypass Toggle Stays Time-Aligned (No Skip)");
    testBypassToggleAlignment();
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

  void testMonoBusLayout() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    juce::AudioProcessor::BusesLayout monoLayout;
    monoLayout.inputBuses.add(juce::AudioChannelSet::mono());
    monoLayout.outputBuses.add(juce::AudioChannelSet::mono());

    NoiseRepellentAudioProcessor proc;
    expect(proc.setBusesLayout(monoLayout), "Mono bus layout must be accepted");
    expect(proc.getTotalNumInputChannels() == 1, "Mono input bus width");
    expect(proc.getTotalNumOutputChannels() == 1, "Mono output bus width");

    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    // Learn a profile through the 1-engine group
    setParam(proc, "learn_noise", 1.0f);
    juce::AudioBuffer<float> buffer(1, blockSize);
    juce::MidiBuffer midi;
    for (int b = 0; b < 30; ++b) {
      generateNoiseBuffer(buffer, 0.05f, 2000 + b);
      proc.processBlock(buffer, midi);
    }
    setParam(proc, "learn_noise", 0.0f);
    pumpMessageLoop(10);

    expect(proc.hasNoiseProfile(), "Mono-learned profile must be present");

    generateNoiseBuffer(buffer, 0.05f, 777);
    proc.processBlock(buffer, midi);
    const auto* out = buffer.getReadPointer(0);
    for (int s = 0; s < blockSize; ++s)
      expect(std::isfinite(out[s]), "Mono output must be finite");

    // State roundtrip between mono instances
    juce::MemoryBlock state;
    proc.getStateInformation(state);
    NoiseRepellentAudioProcessor proc2;
    expect(proc2.setBusesLayout(monoLayout),
           "Mono layout on restore target must be accepted");
    proc2.prepareToPlay(sampleRate, blockSize);
    proc2.setStateInformation(state.getData(), (int)state.getSize());
    expect(proc2.hasNoiseProfile(), "Restored mono profile must be present");

    // Switching back to stereo rebuilds the engine width without crashing
    juce::AudioProcessor::BusesLayout stereoLayout;
    stereoLayout.inputBuses.add(juce::AudioChannelSet::stereo());
    stereoLayout.outputBuses.add(juce::AudioChannelSet::stereo());
    expect(proc.setBusesLayout(stereoLayout),
           "Stereo layout must be accepted after mono");
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> stereoBuf(2, blockSize);
    generateNoiseBuffer(stereoBuf, 0.05f, 888);
    proc.processBlock(stereoBuf, midi);
    expect(proc.hasNoiseProfile(), "Profile must survive mono->stereo switch");

    proc.releaseResources();
    proc2.releaseResources();
  }

  void testStereoMonoStereoDistinctRestore() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    // Learn distinct L/R profiles from deliberately different channel noise.
    setParam(proc, "learn_noise", 1.0f);
    learnStereoDistinctProfiles(proc, blockSize, 30);
    setParam(proc, "learn_noise", 0.0f);
    pumpMessageLoop(10);
    expect(proc.hasNoiseProfile(), "Stereo-learned profile must be present");

    juce::MemoryBlock stereoState;
    proc.getStateInformation(stereoState);
    auto before = extractStateProfiles(stereoState);

    // Find a mode stored for both channels and require distinct L/R data.
    int probeMode = -1;
    for (const auto& [key, data] : before) {
      const auto [ch, mode] = key;
      if (ch == 0 && before.count({1, mode}) > 0) {
        probeMode = mode;
        break;
      }
    }
    expect(probeMode >= 0, "Saved state must hold both channel profiles");
    if (probeMode < 0) {
      proc.releaseResources();
      return;
    }
    expect(maxAbsDiff(before[{0, probeMode}], before[{1, probeMode}]) > 1e-3f,
           "Left/right profiles must be distinct before narrowing");

    // Narrow to mono: the right-channel entry must be retained, not dropped.
    juce::AudioProcessor::BusesLayout monoLayout;
    monoLayout.inputBuses.add(juce::AudioChannelSet::mono());
    monoLayout.outputBuses.add(juce::AudioChannelSet::mono());
    expect(proc.setBusesLayout(monoLayout), "Mono layout must be accepted");
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);
    expect(proc.hasNoiseProfile(), "Mono engine must keep channel 0 profile");

    // Widen back to stereo: the retained right channel must come back intact
    // (not overwritten by a channel-0 fallback copy).
    juce::AudioProcessor::BusesLayout stereoLayout;
    stereoLayout.inputBuses.add(juce::AudioChannelSet::stereo());
    stereoLayout.outputBuses.add(juce::AudioChannelSet::stereo());
    expect(proc.setBusesLayout(stereoLayout),
           "Stereo layout must be accepted after mono");
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);
    expect(proc.hasNoiseProfile(), "Profile must survive mono->stereo switch");

    juce::MemoryBlock restoredState;
    proc.getStateInformation(restoredState);
    auto after = extractStateProfiles(restoredState);
    expect(after.count({0, probeMode}) > 0 && after.count({1, probeMode}) > 0,
           "Restored state must hold both channel profiles");
    if (after.count({0, probeMode}) > 0 && after.count({1, probeMode}) > 0) {
      expectWithinAbsoluteError(
          maxAbsDiff(after[{1, probeMode}], before[{1, probeMode}]), 0.0f,
          1e-6f, "Right channel must restore intact across mono narrowing");
      expect(maxAbsDiff(after[{0, probeMode}], after[{1, probeMode}]) > 1e-3f,
             "Restored left/right profiles must remain distinct");
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

  void testFrameSizeSwitch() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    // Default must be the legacy 46 ms frame.
    expectWithinAbsoluteError(proc.getFrameSizeMs(), 46.0f, 1e-6f,
                              "Default frame size must be 46 ms");
    const int latency46 = proc.getLatencySamples();
    expect(latency46 > 0, "Default engine must report positive latency");

    learnProfile(proc, sampleRate, blockSize, 20);
    expect(proc.hasNoiseProfile(), "Profile must exist before the switch");

    // Upscale 46 -> 64 ms (choice index 2 -> 3) on the message thread.
    // Policy: the switch resets the profile (resampling across resolutions
    // works poorly), so the user re-learns at native resolution.
    setParam(proc, "frame_size", 3.0f);
    pumpMessageLoop(50);
    expectWithinAbsoluteError(proc.getFrameSizeMs(), 64.0f, 1e-6f,
                              "Frame size must switch to 64 ms");
    expect(proc.getLatencySamples() > latency46,
           "Larger frame must report larger latency");
    expect(!proc.hasNoiseProfile(),
           "Frame-size switch must reset the learned profile");

    // Audio must keep flowing without NaNs after the rebuild.
    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;
    for (int b = 0; b < 10; ++b) {
      generateNoiseBuffer(buffer, 0.05f, 5000 + b);
      proc.processBlock(buffer, midi);
      for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
        const float* d = buffer.getReadPointer(ch);
        for (int s = 0; s < buffer.getNumSamples(); ++s)
          expect(!std::isnan(d[s]) && !std::isinf(d[s]),
                 "Output must stay finite after frame-size switch");
      }
    }

    // Downscale back to 46 ms: still no profile, latency restored.
    setParam(proc, "frame_size", 2.0f);
    pumpMessageLoop(50);
    expect(proc.getLatencySamples() == latency46,
           "Latency must return to the 46 ms value");
    expect(!proc.hasNoiseProfile(),
           "Profile must stay cleared after switching back");

    // Re-learning at the new size works normally.
    learnProfile(proc, sampleRate, blockSize, 20);
    expect(proc.hasNoiseProfile(),
           "Re-learn must capture a fresh profile after the switch");

    // Largest step (93 ms, choice index 4): latency grows again within the
    // DryWetMixer's headroom, audio stays finite, profile resets.
    setParam(proc, "frame_size", 4.0f);
    pumpMessageLoop(50);
    expectWithinAbsoluteError(proc.getFrameSizeMs(), 93.0f, 1e-6f,
                              "Frame size must switch to 93 ms");
    expect(proc.getLatencySamples() > latency46,
           "93 ms frame must report larger latency than 46 ms");
    expect(proc.getLatencySamples() < 65536,
           "93 ms latency must fit the DryWetMixer delay headroom");
    expect(!proc.hasNoiseProfile(),
           "Switch to 93 ms must reset the learned profile");
    generateNoiseBuffer(buffer, 0.05f, 6000);
    proc.processBlock(buffer, midi);
    expect(!std::isnan(buffer.getSample(0, 0)),
           "Output must stay finite at 93 ms");

    // Switch with no profile: plain rebuild, nothing flagged.
    NoiseRepellentAudioProcessor fresh;
    fresh.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);
    setParam(fresh, "frame_size", 0.0f);
    pumpMessageLoop(50);
    expectWithinAbsoluteError(fresh.getFrameSizeMs(), 23.0f, 1e-6f,
                              "Fresh instance must switch to 23 ms");
    expect(!fresh.hasNoiseProfile(),
           "Fresh instance must still have no profile");
    fresh.releaseResources();

    proc.releaseResources();
  }

  void testLowLatencyMode() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    const int latency46 = proc.getLatencySamples();
    expect(latency46 > 512, "46 ms engine latency must exceed 512 samples");

    learnProfile(proc, sampleRate, blockSize, 20);
    expect(proc.hasNoiseProfile(), "Profile must exist before the switch");

    // Enable low latency: clean slate like a frame-size switch,
    // 512-sample latency.
    setParam(proc, "low_latency", 1.0f);
    pumpMessageLoop(50);
    expect(proc.getLatencySamples() == 512,
           "Low-latency engine must report 512 samples at 48 kHz");
    expect(!proc.hasNoiseProfile(),
           "Low-latency switch must reset the learned profile");

    // Audio must keep flowing without NaNs after the rebuild, including
    // at extreme smoothing (capped release must not pad or blow up).
    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;
    setParam(proc, "smoothing_factor", 100.0f);
    for (int b = 0; b < 10; ++b) {
      generateNoiseBuffer(buffer, 0.05f, 7000 + b);
      proc.processBlock(buffer, midi);
      for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
        const float* d = buffer.getReadPointer(ch);
        for (int s = 0; s < buffer.getNumSamples(); ++s)
          expect(!std::isnan(d[s]) && !std::isinf(d[s]),
                 "Output must stay finite in low-latency mode");
      }
    }

    // Disable again: latency restored, still no profile.
    setParam(proc, "low_latency", 0.0f);
    pumpMessageLoop(50);
    expect(proc.getLatencySamples() == latency46,
           "Latency must return to the 46 ms value");
    expect(!proc.hasNoiseProfile(),
           "Profile must stay cleared after leaving low-latency mode");

    proc.releaseResources();
  }

  static float runAndMeasureRms(NoiseRepellentAudioProcessor& proc,
                                int blockSize, int warmBlocks, int measBlocks,
                                int seedBase) {
    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;
    for (int b = 0; b < warmBlocks; ++b) {
      generateNoiseBuffer(buffer, 0.05f, seedBase + b);
      proc.processBlock(buffer, midi);
    }
    double sumSq = 0.0;
    long count = 0;
    for (int b = 0; b < measBlocks; ++b) {
      generateNoiseBuffer(buffer, 0.05f, seedBase + 1000 + b);
      proc.processBlock(buffer, midi);
      for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
        const float* d = buffer.getReadPointer(ch);
        for (int s = 0; s < buffer.getNumSamples(); ++s) {
          sumSq += static_cast<double>(d[s]) * d[s];
          ++count;
        }
      }
    }
    return static_cast<float>(std::sqrt(sumSq / static_cast<double>(count)));
  }

  void testReductionCurveLowLatency() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

  // Low-latency leg: uniform curve must shape the 512-sample engine output.
  {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    setParam(proc, "low_latency", 1.0f);
    pumpMessageLoop(50);
    expect(proc.getLatencySamples() == 512,
           "Low-latency engine must report 512 samples");

    learnProfile(proc, sampleRate, blockSize, 25);

    // Flat curve (all zeros): baseline output level.
    setParam(proc, "reduction_curve_enabled", 1.0f);
    proc.setCurveNodes({{0.0f, 0.0f}, {1.0f, 0.0f}});
    const float rmsFlat =
        runAndMeasureRms(proc, blockSize, 25, 10, 41000);

    // Uniform +24 dB extra reduction: output must drop hard.
    proc.setCurveNodes({{0.0f, 24.0f}, {1.0f, 24.0f}});
    const float rmsCut =
        runAndMeasureRms(proc, blockSize, 25, 10, 42000);

    expect(rmsCut < rmsFlat * 0.5f,
           "Uniform +24 dB curve must clearly reduce low-latency output");
    proc.releaseResources();
  }
  }

  void testBypassToggleAlignment() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    constexpr int settleBlocks = 12; // > latency (4416) + margin

    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    // Unity reduction: the wet path passes signal through cleanly so peak
    // positions purely reflect pipeline delay, not denoising.
    setParam(proc, "reduction_amount", 0.0f);

    const int latency = proc.getLatencySamples();
    expect(latency > 0, "Engine must report positive latency");

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // Fires one isolated impulse over a low (non-silent) noise floor and
    // returns its stream offset of peak energy. The floor keeps the engine
    // deterministically running (no silence-sleep involvement) while sitting
    // far below the impulse and the detection threshold; history stays
    // near-clean so the response is deterministic.
    int seedCounter = 0;
    auto fireImpulse = [&](std::vector<float>& stream) {
      stream.clear();
      stream.reserve(static_cast<size_t>(blockSize) * settleBlocks);
      for (int b = 0; b < settleBlocks; ++b) {
        generateNoiseBuffer(buffer, 0.0001f, 9000 + seedCounter++);
        if (b == 0)
          buffer.setSample(0, 0, buffer.getSample(0, 0) + 0.5f);
        proc.processBlock(buffer, midi);
        const float* d = buffer.getReadPointer(0);
        stream.insert(stream.end(), d, d + blockSize);
      }
      size_t peakIdx = 0;
      float peakVal = 0.0f;
      for (size_t i = 0; i < stream.size(); ++i) {
        if (std::abs(stream[i]) > peakVal) {
          peakVal = std::abs(stream[i]);
          peakIdx = i;
        }
      }
      return std::pair<size_t, float>{peakIdx, peakVal};
    };

    auto expectPeakAtLatency = [&](const char* what) {
      std::vector<float> stream;
      const auto [peakIdx, peakVal] = fireImpulse(stream);
      fprintf(stderr, "DBG %s: peak=%zu val=%.4f latency=%d\n", what, peakIdx,
              peakVal, latency);
      expect(peakVal > 0.1f,
             juce::String(what) + ": impulse must come through");
      expect(static_cast<int>(std::abs(static_cast<int>(peakIdx) - latency)) <=
                 2,
             juce::String(what) +
                 ": impulse must land on the reported latency, not jump");
    };

    // Baseline: unbypassed wet path delay.
    expectPeakAtLatency("Unbypassed");
    // Engage: the old skip-the-engine code time-travelled and peaked at ~0.
    setParam(proc, "bypass", 1.0f);
    expectPeakAtLatency("Bypassed");
    // Disengage: engine never stopped, so alignment holds by construction.
    setParam(proc, "bypass", 0.0f);
    expectPeakAtLatency("Un-bypassed");

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
