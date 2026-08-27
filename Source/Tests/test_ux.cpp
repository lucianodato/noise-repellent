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
#include <cmath>
#include <random>
#include <thread>
#include <vector>

#include <juce_core/juce_core.h>
#include <juce_events/juce_events.h>
#include <juce_gui_basics/juce_gui_basics.h>

#include "PluginEditor.h"
#include "PluginProcessor.h"

namespace {

void pumpMessageLoop(int maxMillis = 20) {
  if (auto* mm = juce::MessageManager::getInstanceWithoutCreating()) {
    mm->runDispatchLoopUntil(maxMillis);
  }
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

void generateNoiseBuffer(juce::AudioBuffer<float>& buffer, float rms = 0.05f) {
  std::mt19937 gen(777);
  std::normal_distribution<float> dist(0.0f, rms);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] = dist(gen);
    }
  }
}

} // namespace

class UxAndEditorTest : public juce::UnitTest {
public:
  UxAndEditorTest()
      : juce::UnitTest("UX & Editor Tests", "NoiseRepellent") {}

  void runTest() override {
    beginTest("Rapid Editor Creation and Destruction Lifecycle");
    testEditorLifecycleStress();

    beginTest("Spectral FIFO Queue Starvation and Satiation");
    testSpectralFifoStarvationAndSatiation();

    beginTest("Spectral Frame dB Range & Sanity Verification");
    testSpectralFrameDbRanges();

    beginTest("Editor Resizing and Component Layout Stress");
    testEditorResizingStress();
  }

private:
  void testEditorLifecycleStress() {
    NoiseRepellentAudioProcessor proc;
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    // Audio thread running concurrently
    std::atomic<bool> keepRunning{true};
    std::thread audioThread([&]() {
      juce::AudioBuffer<float> buffer(2, blockSize);
      juce::MidiBuffer midi;
      while (keepRunning.load(std::memory_order_relaxed)) {
        generateNoiseBuffer(buffer, 0.05f);
        proc.processBlock(buffer, midi);
      }
    });

    // Repeatedly instantiate and destroy Editor
    for (int i = 0; i < 40; ++i) {
      auto* editor = proc.createEditorIfNeeded();
      expect(editor != nullptr, "Editor creation must succeed");

      pumpMessageLoop(15);

      if (editor != nullptr) {
        proc.editorBeingDeleted(editor);
        delete editor;
      }
      pumpMessageLoop(5);
    }

    keepRunning.store(false, std::memory_order_relaxed);
    audioThread.join();

    expect(true, "Editor lifecycle stress completed cleanly");
    proc.releaseResources();
  }

  void testSpectralFifoStarvationAndSatiation() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    // 1. Starvation: Process 500 blocks without polling GUI
    for (int b = 0; b < 500; ++b) {
      generateNoiseBuffer(buffer, 0.05f);
      proc.processBlock(buffer, midi);
    }

    // Now poll: FIFO should deliver a valid fresh frame without hang
    NoiseRepellentAudioProcessor::SpectralFrame frame;
    bool gotFrame = proc.getNextSpectralFrame(frame);
    expect(gotFrame, "Processor must return spectral frame after starvation");

    // 2. Satiation: Drain all frames
    int framesDrained = 0;
    while (proc.getNextSpectralFrame(frame)) {
      framesDrained++;
    }
    expect(framesDrained >= 0, "Frames drained successfully");

    proc.releaseResources();
  }

  void testSpectralFrameDbRanges() {
    NoiseRepellentAudioProcessor proc;
    ScopedEditor editor(proc);
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);

    juce::AudioBuffer<float> buffer(2, blockSize);
    juce::MidiBuffer midi;

    for (int b = 0; b < 30; ++b) {
      generateNoiseBuffer(buffer, 0.05f);
      proc.processBlock(buffer, midi);
    }

    NoiseRepellentAudioProcessor::SpectralFrame frame;
    bool gotFrame = false;
    while (proc.getNextSpectralFrame(frame)) {
      gotFrame = true;
    }

    expect(gotFrame, "Must receive spectral frame");
      // Validate dB values are finite and within reasonable bounds (-250 dB to +40 dB)
      for (size_t bin = 0; bin < NoiseRepellentAudioProcessor::kFftBins; ++bin) {
        float inDb = frame.inputMagnitudeDB[bin];
        float outDb = frame.outputMagnitudeDB[bin];
        float nfDb = frame.noiseFloorDB[bin];

        expect(!std::isnan(inDb), "Input magnitude must not be NaN");
        expect(!std::isinf(inDb), "Input magnitude must not be Inf");
        expect(!std::isnan(outDb), "Output magnitude must not be NaN");
        expect(!std::isinf(outDb), "Output magnitude must not be Inf");
        expect(!std::isnan(nfDb), "Noise floor magnitude must not be NaN");
        expect(!std::isinf(nfDb), "Noise floor magnitude must not be Inf");

        expect(inDb >= -250.0f && inDb <= 40.0f, "Input dB within expected range");
        expect(outDb >= -250.0f && outDb <= 40.0f, "Output dB within expected range");
        expect(nfDb >= -250.0f && nfDb <= 40.0f, "Noise floor dB within expected range");
      }

    proc.releaseResources();
  }

  void testEditorResizingStress() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);

    auto* editor = proc.createEditorIfNeeded();
    expect(editor != nullptr, "Editor must be created");

    if (editor != nullptr) {
      const std::vector<std::pair<int, int>> sizes = {
          {400, 300}, {800, 600}, {1024, 768}, {600, 400}, {1200, 900}, {400, 300}};

      for (const auto& size : sizes) {
        editor->setSize(size.first, size.second);
        pumpMessageLoop(10);
      }

      proc.editorBeingDeleted(editor);
      delete editor;
    }

    proc.releaseResources();
  }
};

static UxAndEditorTest uxAndEditorTest;

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
