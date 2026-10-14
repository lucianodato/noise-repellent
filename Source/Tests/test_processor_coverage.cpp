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
#include <cstring>
#include <random>
#include <thread>
#include <vector>

#include <juce_core/juce_core.h>
#include <juce_events/juce_events.h>
#include <juce_gui_basics/juce_gui_basics.h>

#include "PluginProcessor.h"

// Defined in PluginProcessor.cpp; called here to cover the factory path.
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter();

namespace {

void pumpMessageLoop(int maxMillis = 10) {
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

// For worker-thread parameter changes: never pump the message loop off-thread.
void setParamNoPump(NoiseRepellentAudioProcessor& proc,
                    const juce::String& paramId, float value) {
  if (auto* param = proc.getAPVTS().getParameter(paramId)) {
    param->setValueNotifyingHost(param->convertTo0to1(value));
  }
}

void generateNoiseBuffer(juce::AudioBuffer<float>& buffer, float rms = 0.05f,
                         int seed = 42) {
  std::mt19937 gen(static_cast<uint32_t>(seed));
  std::normal_distribution<float> dist(0.0f, rms);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] = dist(gen);
    }
  }
}

void expectAllFinite(juce::UnitTest& t, const juce::AudioBuffer<float>& buffer,
                     const char* msg) {
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    const auto* d = buffer.getReadPointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      if (!std::isfinite(d[s])) {
        t.expect(false, msg);
        return;
      }
    }
  }
  t.expect(true, msg);
}

void processNoise(NoiseRepellentAudioProcessor& proc,
                  juce::AudioBuffer<float>& buffer, int blocks, int seedBase) {
  juce::MidiBuffer midi;
  for (int b = 0; b < blocks; ++b) {
    generateNoiseBuffer(buffer, 0.05f, seedBase + b);
    proc.processBlock(buffer, midi);
  }
}

void learnShortProfile(NoiseRepellentAudioProcessor& proc,
                       int blockSize = 512, int numBlocks = 25,
                       int seedBase = 1000) {
  setParam(proc, "learn_noise", 1.0f);
  juce::AudioBuffer<float> buffer(2, blockSize);
  processNoise(proc, buffer, numBlocks, seedBase);
  setParam(proc, "learn_noise", 0.0f);
}

// Steady two-tone + noise block: drives the tonal-peak detector.
void fillTonalBuffer(juce::AudioBuffer<float>& buffer, int blockIndex) {
  std::mt19937 gen(static_cast<uint32_t>(4242 + blockIndex));
  std::normal_distribution<float> dist(0.0f, 0.008f);
  for (int s = 0; s < buffer.getNumSamples(); ++s) {
    const float t =
        static_cast<float>(blockIndex * buffer.getNumSamples() + s) / 48000.0f;
    const float tone =
        0.5f * std::sin(2.0f * juce::MathConstants<float>::pi * 997.0f * t) +
        0.35f *
            std::sin(2.0f * juce::MathConstants<float>::pi * 3001.0f * t);
    const float n = dist(gen);
    for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
      buffer.setSample(ch, s, tone + n);
    }
  }
}

// Learn a noise profile with strong steady tones so the library reports
// tonal peaks for the profile visualization.
void learnTonalProfile(NoiseRepellentAudioProcessor& proc, int blockSize = 512,
                       int numBlocks = 25) {
  setParam(proc, "learn_noise", 1.0f);
  juce::AudioBuffer<float> buffer(2, blockSize);
  juce::MidiBuffer midi;
  for (int b = 0; b < numBlocks; ++b) {
    fillTonalBuffer(buffer, b);
    proc.processBlock(buffer, midi);
  }
  setParam(proc, "learn_noise", 0.0f);
}

std::vector<float> makeProfileData(size_t n) {
  std::vector<float> d(n);
  for (size_t i = 0; i < n; ++i) {
    d[i] = 0.01f + 0.001f * static_cast<float>(i % 17);
  }
  return d;
}

juce::String base64Of(const std::vector<float>& v) {
  juce::MemoryBlock mb(v.data(), v.size() * sizeof(float));
  return mb.toBase64Encoding();
}

juce::String profileEntryXml(int channel, int mode, int size, int blockCount,
                             const juce::String& dataB64) {
  return "<CHANNEL_PROFILE channel=\"" + juce::String(channel) + "\" mode=\"" +
         juce::String(mode) + "\" size=\"" + juce::String(size) +
         "\" blockCount=\"" + juce::String(blockCount) + "\" data=\"" +
         dataB64 + "\"/>";
}

void restorePlainXml(NoiseRepellentAudioProcessor& proc,
                     const juce::String& xml) {
  proc.setStateInformation(xml.toRawUTF8(),
                           static_cast<int>(xml.getNumBytesAsUTF8()));
}

struct FakePlayHead : public juce::AudioPlayHead {
  bool hasPosition = true;
  bool playing = true;

  juce::Optional<juce::AudioPlayHead::PositionInfo> getPosition()
      const override {
    if (!hasPosition) {
      return {};
    }
    juce::AudioPlayHead::PositionInfo info;
    info.setIsPlaying(playing);
    return info;
  }
};

struct ScopedEditor {
  NoiseRepellentAudioProcessor& proc;
  std::unique_ptr<juce::AudioProcessorEditor> editor;
  explicit ScopedEditor(NoiseRepellentAudioProcessor& p)
      : proc(p), editor(p.createEditorIfNeeded()) {
  }
  ~ScopedEditor() {
    if (editor != nullptr) {
      proc.editorBeingDeleted(editor.get());
    }
  }
};

} // namespace

class ProcessorCoverageTest : public juce::UnitTest {
public:
  ProcessorCoverageTest()
      : juce::UnitTest("Processor Coverage Tests", "NoiseRepellent") {}

  void runTest() override {
    beginTest("Trivial Getters And Unprepared Defaults");
    testTrivialGetters();

    beginTest("Bypass Parameter Per Wrapper Type");
    testBypassParameterWrapperTypes();

    beginTest("Frame Size Getter Across All Choices");
    testFrameSizeGetter();

    beginTest("Parameter Sweep And Combinations");
    testParameterSweep();

    beginTest("Reduction Curve Subsystem Branches");
    testCurveSubsystem();

    beginTest("Curve/Profile Static Decode Fallbacks");
    testStaticDecodeFallbacks();

    beginTest("Profile State Decode Edge Cases");
    testProfileDecodeEdgeCases();

    beginTest("Pending Profile Retain And Rebuild Restore");
    testPendingProfilesRetain();

    beginTest("State Round-Trip: Params, Curve, Profiles");
    testStateRoundTrip();

    beginTest("Corrupt And Plain-XML State Ingestion");
    testCorruptStateIngestion();

    beginTest("Host Transport Variants");
    testTransportVariants();

    beginTest("Non-Realtime Block Counting");
    testNonRealtimeCounting();

    beginTest("Bypass Processing And Zero-Sample Edges");
    testBypassAndZeroSample();

    beginTest("Buses Layout Matrix");
    testBusesLayoutMatrix();

    beginTest("Unprepared Engine And Block-Size Edges");
    testUnpreparedAndBlockSizes();

    beginTest("Silence Streaks And Sleep Gate");
    testSilenceStreaks();

    beginTest("Structural Changes Off The Message Thread");
    testOffThreadStructuralChanges();

    beginTest("Off-Thread Restore Preserves Profile (GH #221 flags)");
    testOffThreadRestorePreserve();

    beginTest("Visualization Feed Direct Paths");
    testVizFeedDirect();

    beginTest("Visualization Through Process Block With Editor");
    testVizThroughProcessBlock();

    beginTest("Engine Param Cache And Forced Load Failure");
    testEngineParamCache();

    beginTest("Editor And Plugin Filter Factory");
    testFactoryAndEditor();
  }

private:
  void testTrivialGetters() {
    NoiseRepellentAudioProcessor proc;
    expectEquals(proc.getName(), juce::String("Noise Repellent"),
                 "getName must report the plugin name");
    expect(!proc.acceptsMidi(), "plugin must not accept MIDI");
    expect(!proc.producesMidi(), "plugin must not produce MIDI");
    expect(!proc.isMidiEffect(), "plugin is not a MIDI effect");
    expectEquals(proc.getTailLengthSeconds(), 0.0, "no tail");
    expectEquals(proc.getNumPrograms(), 1, "single program");
    expectEquals(proc.getCurrentProgram(), 0, "program 0 selected");
    proc.setCurrentProgram(7);
    expectEquals(proc.getCurrentProgram(), 0,
                 "single-state plugin ignores program selection");
    expect(proc.getProgramName(0).isEmpty(), "no program names");
    proc.changeProgramName(0, "renamed");
    expect(proc.hasEditor(), "editor must be available");
    expect(proc.getAPVTS().getParameter("frame_size") != nullptr,
           "APVTS must expose the frame_size parameter");
    expectEquals(proc.getFrameSizeMs(),
                 NoiseRepellentAudioProcessor::kFrameSizeOptionsMs[2],
                 "default frame size is the 46 ms option");
    expect(!proc.isLowLatency(), "low latency off by default");
    expect(!proc.hasNoiseProfile(),
           "unprepared engine must report no noise profile");
    proc.resetNoiseProfile(); // null engine group: must stay safe
    expect(!proc.hasNoiseProfile(), "still no profile after reset");
    expect(proc.engine.group() == nullptr, "engine is not built yet");
    expectEquals(proc.engine.transientIntensity(), 0.0f,
                 "no engine -> zero transient intensity");
    expect(proc.status.isHostTransportPlaying(), "transport defaults to true");
    expect(!proc.status.isHostTransportStateKnown(),
           "transport state unknown before any playhead");
    expectEquals(static_cast<int>(proc.status.getNonRealtimeBlockCount()), 0,
                 "no offline blocks yet");
    expectEquals(proc.status.consumeTransientActivity(), 0.0f,
                 "no transient activity yet");
    expectEquals(proc.status.currentFrameSizeMs, 46.0f,
                 "status mirrors the default frame size");
    expect(proc.profiles.entries().empty(), "no pending profiles");
    proc.profiles.clear();
    {
      juce::SpinLock::ScopedLockType g(proc.profiles.lock());
      juce::ignoreUnused(g);
    }
    {
      juce::SpinLock::ScopedLockType g(proc.curve.lock());
      juce::ignoreUnused(g);
    }
    expectEquals(static_cast<int>(proc.curve.nodes().size()), 2,
                 "curve starts with the two boundary nodes");
    expectEquals(static_cast<int>(proc.curve.biasSize()), 0,
                 "no interpolated bias before any interpolation");

    NoiseRepellentAudioProcessor::SpectralFrame frame;
    expect(!proc.viz.getNextFrame(frame), "empty ring buffer yields no frame");
    proc.viz.resetAccumCount();
    proc.viz.resetAccumulators();
    proc.viz.prepareFrames();
    proc.viz.resizeDelayLine(1024);
  }

  void testBypassParameterWrapperTypes() {
    juce::AudioProcessor::setTypeOfNextNewPlugin(
        juce::AudioProcessor::wrapperType_LV2);
    {
      NoiseRepellentAudioProcessor lv2Proc;
      expect(lv2Proc.getBypassParameter() == nullptr,
             "LV2 wrapper must not receive the bypass parameter");
    }
    juce::AudioProcessor::setTypeOfNextNewPlugin(
        juce::AudioProcessor::wrapperType_Undefined);
    {
      NoiseRepellentAudioProcessor proc;
      expect(proc.getBypassParameter() == proc.getAPVTS().getParameter("bypass"),
             "other wrappers get the internal bypass parameter");
    }
  }

  void testFrameSizeGetter() {
    NoiseRepellentAudioProcessor proc;
    for (int i = 0; i < 5; ++i) {
      setParam(proc, "frame_size", static_cast<float>(i));
      expectEquals(proc.getFrameSizeMs(),
                   NoiseRepellentAudioProcessor::kFrameSizeOptionsMs[
                       static_cast<size_t>(i)],
                   "frame size getter must map every choice index");
    }
    setParam(proc, "low_latency", 1.0f);
    expect(proc.isLowLatency(), "low_latency parameter must be readable");
    setParam(proc, "low_latency", 0.0f);
    expect(!proc.isLowLatency(), "low_latency off must be readable");
  }

  void testParameterSweep() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;

    auto step = [&](const char* id, float value) {
      setParam(proc, id, value);
      generateNoiseBuffer(buffer, 0.05f, 31);
      proc.processBlock(buffer, midi);
      expectAllFinite(*this, buffer, "output must stay finite for parameter");
    };

    // Smoothing method x smoothing amount (smoothingModeFor all branches)
    for (float algo = 0.0f; algo < 4.0f; algo += 1.0f) {
      step("algorithm_mode", algo);
      step("smoothing_factor", 50.0f);
      step("smoothing_factor", 0.0f);
    }
    // Linked vs independent tonal reduction and threshold offsets
    step("link_reduction", 0.0f);
    step("tonal_reduction", 25.0f);
    step("link_threshold_offset", 0.0f);
    step("tonal_noise_profile_offset", -7.0f);
    step("noise_profile_offset", 12.0f);
    step("noise_profile_offset", -12.0f);
    step("link_reduction", 1.0f);
    step("link_threshold_offset", 1.0f);
    // Low latency (forced link + halved smoothing + control defaults push)
    step("low_latency", 1.0f);
    step("low_latency", 0.0f);
    // Adaptive estimation methods
    for (float method = 0.0f; method < 3.0f; method += 1.0f) {
      step("adaptive_method", method);
      step("adaptive_noise", 1.0f);
      step("adaptive_noise", 0.0f);
    }
    // Learn path
    step("learn_noise", 1.0f);
    step("learn_noise", 0.0f);
    // Scalar controls at range ends
    step("aggressiveness", -1.0f);
    step("aggressiveness", 1.0f);
    step("reduction_amount", 40.0f);
    step("reduction_amount", 0.0f);
    step("masking_depth", 0.0f);
    step("masking_depth", 100.0f);
    step("whitening_factor", 100.0f);
    step("whitening_factor", 0.0f);
    step("residual_listen", 1.0f);
    step("residual_listen", 0.0f);
    step("transient_protection_enable", 1.0f);
    step("transient_protection_enable", 0.0f);
    // GUI-only flags still go through buildParams' snapshot consumers
    step("show_advanced", 1.0f);
    step("show_tooltips", 0.0f);
    step("show_hud", 0.0f);
    // Reduction curve enabled with the default nodes
    step("reduction_curve_enabled", 1.0f);
    step("reduction_curve_enabled", 0.0f);
    // Internal bypass on and off (engine keeps running while bypassed)
    step("bypass", 1.0f);
    step("bypass", 0.0f);
    // Frame-size change while low latency is already active: no transition
    // control defaults apply on either side
    step("low_latency", 1.0f);
    step("frame_size", 3.0f);
    step("low_latency", 0.0f);
  }

  void testCurveSubsystem() {
    NoiseRepellentAudioProcessor proc;
    auto& curve = proc.curve;

    // Boundary insertion branches in setNodes
    curve.setNodes({});
    expectEquals(static_cast<int>(curve.nodes().size()), 2,
                 "empty node list gets both boundary nodes");
    curve.setNodes({{0.5f, 5.0f}});
    expectEquals(static_cast<int>(curve.nodes().size()), 3,
                 "interior-only node gets both boundary nodes");
    expectEquals(curve.nodes().front().normX, 0.0f, "front boundary is 0");
    expectEquals(curve.nodes().back().normX, 1.0f, "back boundary is 1");
    curve.setNodes({{0.0f, 1.0f}, {1.0f, -1.0f}});
    expectEquals(static_cast<int>(curve.nodes().size()), 2,
                 "explicit boundaries are kept as-is");

    // addNode clamps and sorts (duplicates exercise zero-width segments)
    curve.addNode(0.5f, 100.0f);
    expectEquals(curve.nodes()[1].normX, 0.5f,
                 "added node sorts between the boundaries");
    expectEquals(curve.nodes()[1].biasDB, 24.0f, "bias clamps to +24 dB");
    curve.addNode(0.25f, -100.0f);
    expectEquals(curve.nodes()[1].normX, 0.25f, "new node keeps its position");
    expectEquals(curve.nodes()[1].biasDB, -24.0f, "bias clamps to -24 dB");
    curve.addNode(0.5f, 1.0f);
    curve.addNode(0.5f, -1.0f);

    // removeNode refuses boundaries and out-of-range indices
    const int before = static_cast<int>(curve.nodes().size());
    curve.removeNode(0);
    curve.removeNode(before - 1);
    curve.removeNode(-3);
    curve.removeNode(before + 5);
    expectEquals(static_cast<int>(curve.nodes().size()), before,
                 "boundary/out-of-range removals are no-ops");
    curve.removeNode(1);
    expectEquals(static_cast<int>(curve.nodes().size()), before - 1,
                 "interior node is removed");

    // updateNode: first/last pinned, middle clamped, out-of-range ignored
    curve.setNodes({{0.0f, 0.0f}, {0.5f, 5.0f}, {1.0f, 0.0f}});
    curve.updateNode(0, 0.9f, 50.0f);
    expectEquals(curve.nodes().front().normX, 0.0f,
                 "first node normX is pinned to 0");
    expectEquals(curve.nodes().front().biasDB, 24.0f, "bias clamps to +24 dB");
    curve.updateNode(2, 0.1f, -50.0f);
    expectEquals(curve.nodes().back().normX, 1.0f,
                 "last node normX is pinned to 1");
    expectEquals(curve.nodes().back().biasDB, -24.0f, "bias clamps to -24 dB");
    curve.updateNode(1, 5.0f, 10.0f);
    expectEquals(curve.nodes()[1].normX, 0.999f,
                 "interior normX clamps below 1");
    curve.updateNode(-1, 0.5f, 1.0f);
    curve.updateNode(99, 0.5f, 1.0f);

    // interpolate edge cases
    curve.interpolate(0, 48000.0); // zero bins: no-op
    curve.interpolate(128, 0.0);  // sr fallback branch
    expectEquals(static_cast<int>(curve.biasSize()), 128,
                 "interpolation fills the requested bin count");
    curve.interpolate(128, -5.0); // sr fallback branch (negative)
    // Multi-node slopes: sign change -> zero slope, same sign -> averaged
    curve.setNodes({{0.0f, 0.0f},
                    {0.3f, 10.0f},
                    {0.6f, -10.0f},
                    {1.0f, 0.0f}});
    curve.interpolate(64, 48000.0);
    curve.setNodes(
        {{0.0f, 0.0f}, {0.3f, 5.0f}, {0.6f, 10.0f}, {1.0f, 20.0f}});
    curve.interpolate(64, 48000.0);
    expectEquals(static_cast<int>(curve.biasSize()), 64, "bin count honored");

    // refreshIfDirty gating
    curve.refreshIfDirty(false, nullptr, 48000.0); // disabled: early return
    curve.setNodes({{0.4f, 2.0f}});
    curve.refreshIfDirty(true, nullptr, 48000.0); // dirty, no group: no resize
    {
      // Lock contention re-arms the dirty flag instead of interpolating
      curve.setNodes({{0.6f, 1.0f}});
      juce::SpinLock::ScopedLockType guard(curve.lock());
      curve.refreshIfDirty(true, nullptr, 48000.0);
    }
    curve.refreshIfDirty(true, nullptr, 48000.0); // consumes the re-armed flag
    curve.refreshIfDirty(true, nullptr, 48000.0); // genuinely clean: early out

    // With a live group the refresh interpolates to the profile size
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    curve.setNodes({{0.5f, 3.0f}});
    curve.refreshIfDirty(true, proc.engine.group(), 48000.0);
    expect(curve.biasSize() > 0, "refresh with engine must interpolate");

    // Spline edge shapes: unsorted input, zero-width segments at the front,
    // at the back, and in the middle (slope guards and t=0 fallback)
    curve.setNodes({{1.0f, 0.0f}, {0.0f, 5.0f}});
    curve.interpolate(64, 48000.0);
    curve.setNodes({{0.0f, 1.0f}, {0.0f, 2.0f}});
    curve.interpolate(64, 48000.0);
    curve.setNodes({{0.0f, 0.0f}, {0.9995f, 1.0f}, {0.9995f, 2.0f}});
    curve.interpolate(64, 48000.0);
    curve.setNodes(
        {{0.0f, 0.0f}, {0.3f, 5.0f}, {0.3f, 9.0f}, {1.0f, 0.0f}});
    curve.interpolate(64, 48000.0);
    expectEquals(static_cast<int>(curve.biasSize()), 64,
                 "edge-shape interpolation must keep the bin count");

    // resetNodes restores the default two-node curve
    curve.resetNodes();
    expectEquals(static_cast<int>(curve.nodes().size()), 2,
                 "resetNodes restores the boundary pair");
    expectEquals(curve.nodes().front().biasDB, 0.0f, "reset bias is flat");
  }

  void testStaticDecodeFallbacks() {
    auto xmlNoCurve = juce::XmlDocument::parse(
        "<PARAMETERS><PARAM id=\"x\" value=\"1\"/></PARAMETERS>");
    auto xmlCurve = juce::XmlDocument::parse(
        "<PARAMETERS><REDUCTION_CURVE><NODE x=\"0.25\" y=\"3\"/><NODE "
        "x=\"0.75\" y=\"-2\"/></REDUCTION_CURVE></PARAMETERS>");
    auto xmlEmptyCurve =
        juce::XmlDocument::parse("<PARAMETERS><REDUCTION_CURVE/></PARAMETERS>");
    auto xmlProfiles = juce::XmlDocument::parse(
        "<PARAMETERS><LEARNED_PROFILES><CHANNEL_PROFILE channel=\"0\" "
        "mode=\"1\" size=\"2\" blockCount=\"1\" data=\"AAAA\"/>"
        "</LEARNED_PROFILES></PARAMETERS>");
    expect(xmlNoCurve != nullptr && xmlCurve != nullptr &&
               xmlEmptyCurve != nullptr && xmlProfiles != nullptr,
           "test XML must parse");

    std::vector<ReductionCurveState::CurveNode> out;

    // Invalid state tree, XML without a curve: both fallbacks exhausted
    expect(!ReductionCurveState::decodeNodesFrom(juce::ValueTree(), *xmlNoCurve,
                                                out),
           "no curve anywhere must decode as absent");
    expect(out.empty(), "nothing decoded");

    // Invalid state tree, XML with a curve: raw-XML fallback
    out.clear();
    expect(ReductionCurveState::decodeNodesFrom(juce::ValueTree(), *xmlCurve,
                                                out),
           "XML fallback must find the curve");
    expectEquals(static_cast<int>(out.size()), 2, "two nodes decoded");
    expectEquals(out[0].normX, 0.25f, "node x parsed");
    expectEquals(out[1].biasDB, -2.0f, "node y parsed");

    // Valid tree without the child, XML with it: fallback again
    const juce::ValueTree stateNoCurve = juce::ValueTree::fromXml(*xmlNoCurve);
    out.clear();
    expect(ReductionCurveState::decodeNodesFrom(stateNoCurve, *xmlCurve, out),
           "fallback must apply when the tree lacks the child");
    expectEquals(static_cast<int>(out.size()), 2, "two nodes decoded");

    // Tree with the child, XML without: tree path
    const juce::ValueTree stateCurve = juce::ValueTree::fromXml(*xmlCurve);
    out.clear();
    expect(ReductionCurveState::decodeNodesFrom(stateCurve, *xmlNoCurve, out),
           "state tree must provide the curve");
    expectEquals(static_cast<int>(out.size()), 2, "two nodes decoded");

    // Empty curve child: rejected
    out.clear();
    expect(!ReductionCurveState::decodeNodesFrom(
               juce::ValueTree::fromXml(*xmlEmptyCurve), *xmlEmptyCurve, out),
           "empty curve child must decode as absent");

    // Node child without x/y attributes: property defaults apply
    auto xmlBareNode = juce::XmlDocument::parse(
        "<PARAMETERS><REDUCTION_CURVE><NODE/></REDUCTION_CURVE></PARAMETERS>");
    expect(xmlBareNode != nullptr, "bare-node XML must parse");
    out.clear();
    expect(ReductionCurveState::decodeNodesFrom(
               juce::ValueTree::fromXml(*xmlBareNode), *xmlBareNode, out),
           "bare node must still decode");
    expectEquals(static_cast<int>(out.size()), 1, "one bare node decoded");
    expectEquals(out[0].normX, 0.0f, "missing x falls back to 0");
    expectEquals(out[0].biasDB, 0.0f, "missing y falls back to 0");

    // findProfilesTree: same fallback matrix
    expect(!ProfileStore::findProfilesTree(juce::ValueTree(), *xmlNoCurve)
                .isValid(),
           "no profiles anywhere must decode as absent");
    expect(ProfileStore::findProfilesTree(juce::ValueTree(), *xmlProfiles)
               .isValid(),
           "XML fallback must find the profiles");
    expect(ProfileStore::findProfilesTree(stateNoCurve, *xmlProfiles)
               .isValid(),
           "fallback must apply when the tree lacks the child");
    expect(ProfileStore::findProfilesTree(
               juce::ValueTree::fromXml(*xmlProfiles), *xmlNoCurve)
               .isValid(),
           "state tree must provide the profiles");
  }

  void testProfileDecodeEdgeCases() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);

    const std::vector<float> data = makeProfileData(32);
    const juce::String goodB64 = base64Of(data);
    const juce::String tinyB64 =
        base64Of(std::vector<float>{1.0f, 2.0f}); // 8 bytes < 32 floats

    auto restoreEntries = [&](const juce::String& entries) {
      restorePlainXml(proc, "<PARAMETERS><LEARNED_PROFILES>" + entries +
                                "</LEARNED_PROFILES></PARAMETERS>");
    };
    auto reset = [&] {
      proc.resetNoiseProfile();
      proc.profiles.clear();
      expect(!proc.hasNoiseProfile(), "reset must clear the profile");
    };

    // Malformed entries must be rejected without touching the engine
    restoreEntries(profileEntryXml(0, 1, 0, 5, goodB64));
    expect(!proc.hasNoiseProfile(), "zero-size entry must be rejected");
    reset();
    restoreEntries(profileEntryXml(0, 1, 32, 5, ""));
    expect(!proc.hasNoiseProfile(), "empty data entry must be rejected");
    reset();
    restoreEntries(profileEntryXml(0, 0, 32, 5, goodB64));
    expect(!proc.hasNoiseProfile(), "mode 0 entry must be rejected");
    reset();
    restoreEntries(profileEntryXml(0, 5, 32, 5, goodB64));
    expect(!proc.hasNoiseProfile(), "out-of-range mode entry must be rejected");
    reset();
    restoreEntries(profileEntryXml(0, 1, 32, 5, tinyB64));
    expect(!proc.hasNoiseProfile(),
           "entry with data shorter than size must be rejected");
    reset();
    restoreEntries(profileEntryXml(0, 1, 32, 5, "!@#$%^&*()"));
    expect(!proc.hasNoiseProfile(),
           "entry with invalid base64 must be rejected");
    reset();

    // Undecodable-for-the-engine entries are retained as pending
    restoreEntries(profileEntryXml(-1, 1, 32, 5, goodB64));
    expect(!proc.hasNoiseProfile(), "negative channel cannot load");
    expect(!proc.profiles.entries().empty(),
           "unloadable entry must be retained as pending");
    reset();
    restoreEntries(profileEntryXml(7, 1, 32, 5, goodB64));
    expect(!proc.hasNoiseProfile(), "channel beyond width cannot load");
    expect(!proc.profiles.entries().empty(),
           "wide entry must be retained for a later rebuild");
    reset();

    // Valid channel-0-only entry: loads, and the sync fills channel 1
    restoreEntries(profileEntryXml(0, 1, 32, 5, goodB64));
    expect(proc.hasNoiseProfile(), "valid entry must load into the engine");
    {
      juce::MemoryBlock state;
      proc.getStateInformation(state);
      std::unique_ptr<juce::XmlElement> xml(
          juce::AudioProcessor::getXmlFromBinary(state.getData(),
                                                 (int)state.getSize()));
      expect(xml != nullptr, "state must serialize to XML");
      if (xml != nullptr) {
        const juce::ValueTree tree = juce::ValueTree::fromXml(*xml);
        // appendToState appends the fresh tree; a LEARNED_PROFILES child
        // carried over from an earlier restore shadows it for findProfilesTree,
        // so assert on the freshly appended (last) tree.
        juce::ValueTree fresh;
        for (int i = 0; i < tree.getNumChildren(); ++i) {
          if (tree.getChild(i).hasType("LEARNED_PROFILES")) {
            fresh = tree.getChild(i);
          }
        }
        expect(fresh.isValid() && fresh.getNumChildren() >= 2,
               "legacy single-channel profile must serialize both channels");
      }
    }

    // Mixed state: junk entry retained while the good one loads
    reset();
    restoreEntries(profileEntryXml(7, 1, 32, 5, goodB64) +
                   profileEntryXml(0, 1, 32, 5, goodB64));
    expect(proc.hasNoiseProfile(), "good entry loads alongside junk");
    expectEquals(static_cast<int>(proc.profiles.entries().size()), 1,
                 "only the junk entry stays pending");

    // captureLive round-trips the live profiles
    std::vector<ProfileStore::SavedProfile> saved;
    ProfileStore::captureLive(proc.engine.group(), saved);
    expect(!saved.empty(), "live profiles must be capturable");
    proc.releaseResources();
  }

  void testPendingProfilesRetain() {
    NoiseRepellentAudioProcessor proc;
    const std::vector<float> data = makeProfileData(32);
    const juce::String goodB64 = base64Of(data);

    // Restore before any engine exists: everything is retained pending
    restorePlainXml(proc, "<PARAMETERS><LEARNED_PROFILES>" +
                              profileEntryXml(7, 1, 32, 5, goodB64) +
                              profileEntryXml(0, 1, 32, 5, goodB64) +
                              "</LEARNED_PROFILES></PARAMETERS>");
    expect(!proc.hasNoiseProfile(), "no engine yet: nothing loaded");
    expectEquals(static_cast<int>(proc.profiles.entries().size()), 2,
                 "both entries retained across the unprepared restore");

    // First engine build takes what fits and keeps the rest
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    expect(proc.hasNoiseProfile(),
           "pending profile must load on the first engine build");
    expectEquals(static_cast<int>(proc.profiles.entries().size()), 1,
                 "the too-wide entry stays pending");

    // Retained entries + no channel-1 data: the legacy fallback loop runs
    // over the pending list during the next restore pass
    std::vector<ProfileStore::SavedProfile> none;
    ProfileStore::PendingProfile bad;
    bad.channel = 0;
    bad.mode = 9; // invalid mode: cannot load, must stay retained
    bad.size = 32;
    bad.blockCount = 5;
    bad.data = data;
    proc.profiles.entries().push_back(bad);
    proc.profiles.restoreIntoNewGroup(proc.engine.group(), none);
    expectEquals(static_cast<int>(proc.profiles.entries().size()), 2,
                 "unloadable entry must survive the rebuild pass");

    // Same legacy fallback through restoreFromTree (ch0-only decoded tree)
    const juce::ValueTree tree = juce::ValueTree::fromXml(
        *juce::XmlDocument::parse("<LEARNED_PROFILES>" +
                                  profileEntryXml(0, 1, 32, 5, goodB64) +
                                  "</LEARNED_PROFILES>"));
    proc.profiles.restoreFromTree(tree, proc.engine.group());
    expect(proc.hasNoiseProfile(), "decoded tree must load");
    expectEquals(static_cast<int>(proc.profiles.entries().size()), 2,
                 "still-pending entries are kept");

    // A later stereo rebuild restores and then clears what fits
    proc.releaseResources();
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    expect(proc.hasNoiseProfile(), "profile must survive re-prepare");

    // decodeProfileEntries over an empty tree must be a clean no-op
    restorePlainXml(proc, "<PARAMETERS><LEARNED_PROFILES/></PARAMETERS>");
    expect(proc.hasNoiseProfile(),
           "empty profiles tree must leave the engine untouched");

    // loadProfileResampledChecked guards: null data and zero size fail
    // without touching the engine and stay retained
    ProfileStore::PendingProfile emptyData;
    emptyData.channel = 0;
    emptyData.mode = 1;
    emptyData.size = 32; // claims data but the vector is empty (null pointer)
    emptyData.blockCount = 5;
    ProfileStore::PendingProfile zeroSize;
    zeroSize.channel = 0;
    zeroSize.mode = 1;
    zeroSize.size = 0;
    zeroSize.data = data;
    proc.profiles.entries().push_back(emptyData);
    proc.profiles.entries().push_back(zeroSize);
    std::vector<ProfileStore::SavedProfile> none2;
    proc.profiles.restoreIntoNewGroup(proc.engine.group(), none2);
    expectEquals(static_cast<int>(proc.profiles.entries().size()), 2,
                 "guarded entries must be retained, not loaded");

    // Serializing a null group produces an empty profiles tree
    juce::ValueTree orphan("PARAMETERS");
    proc.profiles.appendToState(nullptr, orphan);
    const juce::ValueTree orphanProfiles =
        orphan.getChildWithName("LEARNED_PROFILES");
    expect(orphanProfiles.isValid() && orphanProfiles.getNumChildren() == 0,
           "null-group serialization must yield an empty tree");

    // Pending entries with a channel-1 member: the legacy ch0->ch1 fallback
    // is skipped entirely
    ProfileStore::PendingProfile ch1;
    ch1.channel = 1;
    ch1.mode = 9; // invalid so it stays pending
    ch1.size = 32;
    ch1.blockCount = 5;
    ch1.data = data;
    proc.profiles.entries().push_back(ch1);
    proc.profiles.restoreIntoNewGroup(proc.engine.group(), none2);
    proc.profiles.restoreFromTree(
        juce::ValueTree::fromXml(
            *juce::XmlDocument::parse("<LEARNED_PROFILES>" +
                                      profileEntryXml(1, 1, 32, 5, goodB64) +
                                      "</LEARNED_PROFILES>")),
        proc.engine.group());
    expect(proc.hasNoiseProfile(), "ch1-carrying restore must still load");

    // Mono engine group: the cross-channel fallback never applies
    NoiseRepellentAudioProcessor monoProc;
    juce::AudioProcessor::BusesLayout monoLayout;
    monoLayout.inputBuses.add(juce::AudioChannelSet::mono());
    monoLayout.outputBuses.add(juce::AudioChannelSet::mono());
    expect(monoProc.setBusesLayout(monoLayout), "mono layout accepted");
    monoProc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    ProfileStore::PendingProfile monoPending;
    monoPending.channel = 0;
    monoPending.mode = 9;
    monoPending.size = 32;
    monoPending.blockCount = 5;
    monoPending.data = data;
    monoProc.profiles.entries().push_back(monoPending);
    std::vector<ProfileStore::SavedProfile> none3;
    monoProc.profiles.restoreIntoNewGroup(monoProc.engine.group(), none3);
    monoProc.profiles.restoreFromTree(
        juce::ValueTree::fromXml(
            *juce::XmlDocument::parse("<LEARNED_PROFILES>" +
                                      profileEntryXml(0, 1, 32, 5, goodB64) +
                                      "</LEARNED_PROFILES>")),
        monoProc.engine.group());
    expect(monoProc.hasNoiseProfile(), "mono restore must load ch0");
    monoProc.releaseResources();

    // Sample-rate change rebuilds the engine in place
    juce::MidiBuffer midi;
    proc.prepareToPlay(44100.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> at44k(2, 512);
    generateNoiseBuffer(at44k, 0.05f, 47);
    proc.processBlock(at44k, midi);
    expectAllFinite(*this, at44k, "rate-change rebuild must process");

    // Channel-count change rebuilds without a rate change
    expect(proc.setBusesLayout(monoLayout), "mono relayout accepted");
    proc.prepareToPlay(44100.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> monoBuf(1, 512);
    generateNoiseBuffer(monoBuf, 0.05f, 48);
    proc.processBlock(monoBuf, midi);
    expectAllFinite(*this, monoBuf, "relayout rebuild must process");
    juce::AudioProcessor::BusesLayout stereoLayout;
    stereoLayout.inputBuses.add(juce::AudioChannelSet::stereo());
    stereoLayout.outputBuses.add(juce::AudioChannelSet::stereo());
    expect(proc.setBusesLayout(stereoLayout), "stereo relayout accepted");
    proc.prepareToPlay(44100.0, 512);
    pumpMessageLoop(10);

    // Explicit default construction of the snapshot types
    ProfileStore::SavedProfile savedDefault;
    ProfileStore::PendingProfile pendingDefault;
    expectEquals(savedDefault.channel, 0, "snapshot defaults");
    expectEquals(pendingDefault.mode, 1, "pending defaults");
    proc.releaseResources();
  }

  void testStateRoundTrip() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    NoiseRepellentAudioProcessor procA;
    procA.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);
    learnShortProfile(procA);
    expect(procA.hasNoiseProfile(), "profile must be learned");
    procA.curve.setNodes(
        {{0.2f, 3.0f}, {0.5f, -1.5f}, {0.8f, -2.0f}});
    setParam(procA, "aggressiveness", 0.8f);
    setParam(procA, "algorithm_mode", 1.0f);
    setParam(procA, "smoothing_factor", 40.0f);
    setParam(procA, "link_threshold_offset", 0.0f);
    setParam(procA, "tonal_noise_profile_offset", 2.5f);

    juce::MemoryBlock state;
    procA.getStateInformation(state);
    expect(state.getSize() > 0, "state must serialize");

    NoiseRepellentAudioProcessor procB;
    procB.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);
    procB.setStateInformation(state.getData(), (int)state.getSize());
    pumpMessageLoop(10);

    expect(procB.hasNoiseProfile(), "restored state must carry the profile");
    expectEquals(static_cast<int>(procB.curve.nodes().size()),
                 static_cast<int>(procA.curve.nodes().size()),
                 "curve node count must round-trip");
    for (size_t i = 0; i < procA.curve.nodes().size(); ++i) {
      expectEquals(procB.curve.nodes()[i].normX,
                   procA.curve.nodes()[i].normX, "node x round-trips");
      expectEquals(procB.curve.nodes()[i].biasDB,
                   procA.curve.nodes()[i].biasDB, "node y round-trips");
    }
    expect(std::abs(procB.getAPVTS()
                        .getRawParameterValue("aggressiveness")
                        ->load() -
                    0.8f) < 0.05f,
           "parameter value must round-trip");
    expect(std::abs(procB.getAPVTS()
                        .getRawParameterValue("tonal_noise_profile_offset")
                        ->load() -
                    2.5f) < 0.05f,
           "unlinked tonal offset must round-trip");
    procA.releaseResources();
    procB.releaseResources();
  }

  void testCorruptStateIngestion() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);

    // Null state
    proc.setStateInformation(nullptr, 0);

    // Non-XML binary garbage: the plain-XML fallback must refuse it
    std::vector<char> garbage(256);
    for (size_t i = 0; i < garbage.size(); ++i) {
      garbage[i] = static_cast<char>(0x80 + (i % 64));
    }
    proc.setStateInformation(garbage.data(), (int)garbage.size());
    // Non-null data with a zero size must also be refused
    proc.setStateInformation(garbage.data(), 0);

    // Truncated XML
    const char* truncXml = "<PARAMETERS><VALUE id=\"aggressiveness\"";
    proc.setStateInformation(truncXml, (int)std::strlen(truncXml));

    // Plain XML text (not the binary envelope) with parameters and a curve
    const juce::String plain =
        "<PARAMETERS><PARAM id=\"aggressiveness\" value=\"0.7\"/>"
        "<PARAM id=\"reduction_amount\" value=\"20\"/>"
        "<REDUCTION_CURVE><NODE x=\"0.1\" y=\"2\"/><NODE x=\"0.9\" y=\"-3\"/>"
        "</REDUCTION_CURVE></PARAMETERS>";
    restorePlainXml(proc, plain);
    expect(std::abs(proc.getAPVTS()
                        .getRawParameterValue("aggressiveness")
                        ->load() -
                    0.7f) < 0.05f,
           "plain-XML state must restore parameter values");
    expectEquals(static_cast<int>(proc.curve.nodes().size()), 4,
                 "plain-XML state must restore the curve nodes");

    // Plain XML carrying a valid profile
    restorePlainXml(proc, "<PARAMETERS><LEARNED_PROFILES>" +
                              profileEntryXml(0, 1, 32, 5,
                                              base64Of(makeProfileData(32))) +
                              "</LEARNED_PROFILES></PARAMETERS>");
    expect(proc.hasNoiseProfile(), "plain-XML state must restore profiles");

    // Binary state produced by the processor itself
    juce::MemoryBlock state;
    proc.getStateInformation(state);
    proc.resetNoiseProfile();
    proc.setStateInformation(state.getData(), (int)state.getSize());
    expect(proc.hasNoiseProfile(), "binary state must restore profiles");
    proc.releaseResources();
  }

  void testTransportVariants() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> buffer(2, 512);
    generateNoiseBuffer(buffer, 0.05f, 17);
    juce::MidiBuffer midi;

    FakePlayHead stopped;
    stopped.playing = false;
    proc.setPlayHead(&stopped);
    proc.processBlock(buffer, midi);
    expect(proc.status.isHostTransportStateKnown(),
           "stopped host transport must be known");
    expect(!proc.status.isHostTransportPlaying(),
           "stopped host transport must report not playing");

    FakePlayHead playing;
    proc.setPlayHead(&playing);
    proc.processBlock(buffer, midi);
    expect(proc.status.isHostTransportPlaying(),
           "playing host transport must report playing");

    FakePlayHead positionLess;
    positionLess.hasPosition = false;
    proc.setPlayHead(&positionLess);
    proc.processBlock(buffer, midi);
    expect(!proc.status.isHostTransportStateKnown(),
           "playhead without position leaves transport unknown");

    proc.setPlayHead(nullptr);
    proc.processBlock(buffer, midi);
    expect(!proc.status.isHostTransportStateKnown(),
           "no playhead leaves transport unknown");
    proc.releaseResources();
  }

  void testNonRealtimeCounting() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> buffer(2, 512);
    generateNoiseBuffer(buffer, 0.05f, 21);
    juce::MidiBuffer midi;

    proc.setNonRealtime(true);
    proc.processBlock(buffer, midi);
    proc.processBlock(buffer, midi);
    proc.processBlockBypassed(buffer, midi);
    expectEquals(static_cast<int>(proc.status.getNonRealtimeBlockCount()), 3,
                 "offline blocks must be counted across both entry points");
    proc.setNonRealtime(false);
    proc.processBlock(buffer, midi);
    expectEquals(static_cast<int>(proc.status.getNonRealtimeBlockCount()), 3,
                 "realtime blocks must not advance the counter");
    proc.releaseResources();
  }

  void testBypassAndZeroSample() {
    NoiseRepellentAudioProcessor proc;
    juce::MidiBuffer midi;
    juce::AudioBuffer<float> empty(2, 0);
    proc.processBlock(empty, midi); // zero-sample early return
    proc.processBlockBypassed(empty, midi);

    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> buffer(2, 512);

    generateNoiseBuffer(buffer, 0.05f, 33);
    proc.processBlockBypassed(buffer, midi);
    expectAllFinite(*this, buffer, "bypassed processing must stay finite");

    setParam(proc, "bypass", 1.0f);
    generateNoiseBuffer(buffer, 0.05f, 34);
    proc.processBlock(buffer, midi);
    expectAllFinite(*this, buffer, "internal bypass must stay finite");
    setParam(proc, "bypass", 0.0f);
    generateNoiseBuffer(buffer, 0.05f, 35);
    proc.processBlock(buffer, midi);
    expectAllFinite(*this, buffer, "processing must stay finite");
    proc.releaseResources();
  }

  void testBusesLayoutMatrix() {
    NoiseRepellentAudioProcessor proc;
    auto make = [](juce::AudioChannelSet in, juce::AudioChannelSet out) {
      juce::AudioProcessor::BusesLayout layout;
      layout.inputBuses.add(in);
      layout.outputBuses.add(out);
      return layout;
    };
    expect(!proc.isBusesLayoutSupported(
               make(juce::AudioChannelSet::quadraphonic(),
                    juce::AudioChannelSet::quadraphonic())),
           "layouts wider than stereo must be rejected");
    expect(!proc.isBusesLayoutSupported(
               make(juce::AudioChannelSet::mono(),
                    juce::AudioChannelSet::stereo())),
           "input/output width mismatch must be rejected");
    expect(proc.isBusesLayoutSupported(
               make(juce::AudioChannelSet::stereo(),
                    juce::AudioChannelSet::stereo())),
           "stereo must be supported");
    expect(proc.isBusesLayoutSupported(make(juce::AudioChannelSet::mono(),
                                            juce::AudioChannelSet::mono())),
           "mono must be supported");
  }

  void testUnpreparedAndBlockSizes() {
    NoiseRepellentAudioProcessor proc;
    juce::AudioBuffer<float> buffer(2, 512);
    generateNoiseBuffer(buffer, 0.05f, 41);
    juce::MidiBuffer midi;

    // processBlock before any prepareToPlay: engine and mixer stay inert
    // (JUCE debug asserts from the unprepared DryWetMixer are logged, not
    // fatal, and its fifos reject the writes)
    {
      NoiseRepellentAudioProcessor fresh;
      generateNoiseBuffer(buffer, 0.05f, 40);
      fresh.processBlock(buffer, midi);
      fresh.processBlockBypassed(buffer, midi);
    }

    // Engine run without any group: must return without touching memory
    generateNoiseBuffer(buffer, 0.05f, 41);
    const auto ep = proc.engine.buildParams(proc.getAPVTS());
    proc.engine.run(buffer, ep);
    expectAllFinite(*this, buffer, "unprepared engine run must be safe");

    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);

    // Zero-channel buffers hit both early-return guards
    juce::AudioBuffer<float> zeroCh(0, 512);
    proc.processBlock(zeroCh, midi);
    proc.processBlockBypassed(zeroCh, midi);
    expect(true, "zero-channel buffers must early-return");

    // Block smaller than the prepared size
    juce::AudioBuffer<float> small(2, 64);
    generateNoiseBuffer(small, 0.05f, 42);
    proc.processBlock(small, midi);
    expectAllFinite(*this, small, "small block must process");

    // Block larger than the advertised host block size
    juce::AudioBuffer<float> big(2, 4096);
    generateNoiseBuffer(big, 0.05f, 43);
    proc.processBlock(big, midi);
    expectAllFinite(*this, big, "large block must process");

    // Direct engine run with a narrower buffer than the group: the run
    // clamps to the buffer width (no mixer involved)
    juce::AudioBuffer<float> narrow(1, 512);
    generateNoiseBuffer(narrow, 0.05f, 44);
    proc.engine.run(narrow, proc.engine.buildParams(proc.getAPVTS()));
    expectAllFinite(*this, narrow, "narrow engine run must stay finite");

    // Full mono-bus pipeline end to end
    NoiseRepellentAudioProcessor monoProc;
    juce::AudioProcessor::BusesLayout monoLayout;
    monoLayout.inputBuses.add(juce::AudioChannelSet::mono());
    monoLayout.outputBuses.add(juce::AudioChannelSet::mono());
    expect(monoProc.setBusesLayout(monoLayout), "mono layout must be accepted");
    monoProc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> mono(1, 512);
    generateNoiseBuffer(mono, 0.05f, 46);
    monoProc.processBlock(mono, midi);
    expectAllFinite(*this, mono, "mono bus must process");
    monoProc.releaseResources();

    // A frame size below one sample at this rate fails engine construction:
    // the rebuild guard and latency reporting must tolerate a null group
    // with a positive sample rate
    NoiseRepellentAudioProcessor tiny;
    tiny.prepareToPlay(1.0, 512);
    expect(tiny.engine.group() == nullptr,
           "degenerate sample rate must fail the engine build");
    setParam(tiny, "frame_size", 4.0f);
    expectEquals(tiny.getFrameSizeMs(), 93.0f,
                 "parameter must still update without an engine");
    tiny.releaseResources();

    // releaseResources then re-prepare
    proc.releaseResources();
    proc.prepareToPlay(48000.0, 256);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> again(2, 512);
    generateNoiseBuffer(again, 0.05f, 45);
    proc.processBlock(again, midi);
    expectAllFinite(*this, again, "processing must resume after re-prepare");
    proc.releaseResources();
  }

  void testSilenceStreaks() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    juce::MidiBuffer midi;
    juce::AudioBuffer<float> silent(2, 512);
    silent.clear();

    // No profile: after the flush margin the engine may sleep. The block is
    // re-filled every iteration (processBlock overwrites it with output) with
    // 5e-6: silent per the -100 dB gate but non-zero, so the sleep path's
    // buffer.clear() is observable as exact zeros.
    for (int b = 0; b < 200; ++b) {
      silent.clear();
      for (int ch = 0; ch < silent.getNumChannels(); ++ch) {
        silent.setSample(ch, 0, 5e-6f);
      }
      proc.processBlock(silent, midi);
    }
    expect(proc.status.silentBlocksStreak >= 100,
           "silence streak must accumulate across silent blocks");
    {
      const auto* d = silent.getReadPointer(0);
      bool allZero = true;
      for (int s = 0; s < silent.getNumSamples(); ++s) {
        allZero = allZero && d[s] == 0.0f;
      }
      expect(allZero,
             "slept engine must clear the block to exact silence");
    }

    // A non-silent block resets the streak
    generateNoiseBuffer(silent, 0.05f, 51);
    proc.processBlock(silent, midi);
    expectEquals(static_cast<int>(proc.status.silentBlocksStreak), 0,
                 "non-silent input must reset the streak");
    silent.clear();

    // Bypass keeps the engine running even when the sleep gate would fire:
    // the block must keep flowing instead of being cleared to silence
    setParam(proc, "bypass", 1.0f);
    for (int b = 0; b < 200; ++b) {
      silent.clear();
      for (int ch = 0; ch < silent.getNumChannels(); ++ch) {
        silent.setSample(ch, 0, 5e-6f);
      }
      proc.processBlock(silent, midi);
    }
    expectAllFinite(*this, silent, "bypassed silence must stay finite");
    {
      const auto* d = silent.getReadPointer(0);
      bool anyNonZero = false;
      for (int s = 0; s < silent.getNumSamples(); ++s) {
        anyNonZero = anyNonZero || d[s] != 0.0f;
      }
      expect(anyNonZero,
             "bypassed input must flow through instead of being cleared");
    }
    setParam(proc, "bypass", 0.0f);

    // Learning keeps the engine awake on silence (streak never advances)
    setParam(proc, "learn_noise", 1.0f);
    for (int b = 0; b < 12; ++b) {
      silent.clear();
      proc.processBlock(silent, midi);
    }
    expectEquals(static_cast<int>(proc.status.silentBlocksStreak), 0,
                 "learning must not count as silence");
    setParam(proc, "learn_noise", 0.0f);

    // Adaptive tracking keeps the engine awake on silence as well
    setParam(proc, "adaptive_noise", 1.0f);
    for (int b = 0; b < 12; ++b) {
      silent.clear();
      proc.processBlock(silent, midi);
    }
    expectEquals(static_cast<int>(proc.status.silentBlocksStreak), 0,
                 "adaptive tracking must not count as silence");
    setParam(proc, "adaptive_noise", 0.0f);

    // With a learned profile the engine keeps running through silence
    learnShortProfile(proc);
    expect(proc.hasNoiseProfile(), "profile learned before the silent run");
    for (int b = 0; b < 80; ++b) {
      silent.clear();
      proc.processBlock(silent, midi);
    }
    expect(proc.hasNoiseProfile(),
           "silence must not discard the learned profile");
    expectAllFinite(*this, silent, "silent tail must stay finite");
    proc.releaseResources();
  }

  void testOffThreadStructuralChanges() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    learnShortProfile(proc, 512, 12, 6000);
    expect(proc.hasNoiseProfile(), "profile present before the switch");

    // Off-thread automation coalesces into one deferred rebuild
    std::thread automator(
        [&proc] { setParamNoPump(proc, "frame_size", 4.0f); });
    automator.join();
    expect(proc.status.frameSizeRebuildPending.load(),
           "off-thread structural change must arm the deferred rebuild");
    expect(!proc.status.preserveProfilesOnFrameRebuild.load(),
           "plain automation must not arm the profile-preserve gate");

    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    generateNoiseBuffer(buffer, 0.05f, 71);
    proc.processBlock(buffer, midi);
    expect(!proc.status.frameSizeRebuildPending.load(),
           "processBlock must consume the rebuild request");
    expectEquals(proc.getFrameSizeMs(), 93.0f,
                 "deferred rebuild must apply the new frame size");
    expect(!proc.hasNoiseProfile(),
                 "frame-size switch starts clean (profile dropped)");
    proc.releaseResources();
  }

  void testOffThreadRestorePreserve() {
    constexpr double sampleRate = 48000.0;
    constexpr int blockSize = 512;

    juce::MemoryBlock saved;
    {
      NoiseRepellentAudioProcessor donor;
      donor.prepareToPlay(sampleRate, blockSize);
      pumpMessageLoop(10);
      learnShortProfile(donor, blockSize, 25, 8000);
      expect(donor.hasNoiseProfile(), "donor must hold a profile");
      donor.getStateInformation(saved);
      donor.releaseResources();
    }

    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(sampleRate, blockSize);
    pumpMessageLoop(10);
    // Make the live frame size differ from the restored one
    setParam(proc, "frame_size", 0.0f);
    pumpMessageLoop(10);

    std::thread restorer([&proc, &saved] {
      proc.setStateInformation(saved.getData(), (int)saved.getSize());
    });
    restorer.join();

    expect(proc.status.preserveProfilesOnFrameRebuild.load(),
           "off-thread restore must arm the profile-preserve gate [GH #221]");
    expect(proc.status.frameSizeRebuildPending.load(),
           "structural restore change defers the rebuild");
    expect(proc.hasNoiseProfile(), "restored profile present before rebuild");

    juce::AudioBuffer<float> buffer(2, blockSize);
    processNoise(proc, buffer, 20, 9000);
    pumpMessageLoop(20);
    expect(!proc.status.preserveProfilesOnFrameRebuild.load(),
           "the preserve gate must be consumed by the rebuild");
    expect(!proc.status.frameSizeRebuildPending.load(),
           "the deferred rebuild must have run");
    expect(proc.hasNoiseProfile(),
           "restored profile must survive the deferred rebuild [GH #221]");
    proc.releaseResources();
  }

  void testVizFeedDirect() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);

    std::vector<float> input(4096, 0.1f);
    std::vector<float> output(4096, 0.05f);

    // Null engine group: frames still push, the profile query bails out
    VisualizationFeed::FeedInput in{
        .group = nullptr,
        .profileLock = proc.profiles.lock(),
        .apvts = proc.getAPVTS(),
        .latency = 0,
        .isLearning = false,
        .adaptiveNoise = false,
        .transientIntensity = 0.5f,
        .transientProtectionActive = true};
    proc.viz.feed(in, std::span<const float>(input.data(), input.size()), std::span<const float>(output.data(), output.size()));
    NoiseRepellentAudioProcessor::SpectralFrame frame;
    auto drainAll = [&frame, &proc] {
      while (proc.viz.getNextFrame(frame)) {
      }
    };
    expect(proc.viz.getNextFrame(frame), "a full window must push a frame");
    expectEquals(frame.transientIntensity, 0.5f,
                 "transient max in window must be reported");
    expect(frame.isTransientProtected, "0.5 intensity counts as protected");
    expect(frame.isTransientProtectionActive,
           "protection-active flag mirrors the parameter snapshot");
    expect(!frame.hasNoiseProfile, "null group cannot hold a profile");

    // Profile-lock contention: the previous frame's floor is retained
    {
      juce::SpinLock::ScopedLockType guard(proc.profiles.lock());
      proc.viz.feed(in, std::span<const float>(input.data(), 1024), std::span<const float>());
    }
    expect(proc.viz.getNextFrame(frame),
           "contention must still push the spectral frame");
    expect(!frame.hasNoiseProfile,
           "busy profile lock must leave the frame without a profile");

    // Silent frames without any profile
    for (int i = 0; i < 4; ++i) {
      proc.viz.pushSilentFrame(in, proc.profiles);
    }
    expect(proc.viz.getNextFrame(frame),
           "every fourth silent block pushes a frame");
    expect(!frame.hasNoiseProfile, "no profile learned yet");
    expectEquals(frame.noiseFloorDB[0], -120.0f,
                 "missing profile floor falls to -120 dB");
    expectEquals(frame.transientIntensity, 0.0f,
                 "silent frames carry no transient energy");

    // Learn a tonal profile, then cover the frozen-profile silent paths
    learnTonalProfile(proc);
    expect(proc.hasNoiseProfile(), "profile learned for the silent paths");
    // Process tonal material against the learned profile so the engine's
    // tonal detector reports peaks
    {
      juce::AudioBuffer<float> tbuf(2, 512);
      juce::MidiBuffer tmidi;
      for (int b = 0; b < 10; ++b) {
        fillTonalBuffer(tbuf, 500 + b);
        proc.processBlock(tbuf, tmidi);
      }
    }
    VisualizationFeed::FeedInput inGroup{
        .group = proc.engine.group(),
        .profileLock = proc.profiles.lock(),
        .apvts = proc.getAPVTS(),
        .latency = proc.status.currentLatency,
        .isLearning = false,
        .adaptiveNoise = false,
        .transientIntensity = 0.0f,
        .transientProtectionActive = false};

    // Learning-path visualization: the rolling-mean profile is queried
    drainAll();
    VisualizationFeed::FeedInput inLearning = inGroup;
    inLearning.isLearning = true;
    proc.viz.feed(inLearning, std::span<const float>(input.data(), input.size()), std::span<const float>(output.data(), output.size()));
    drainAll();
    expect(frame.hasNoiseProfile, "learning frame sees the rolling profile");
    expect(frame.hasNoiseFloorEstimate,
           "learning frame carries the floor estimate");

    // Full FFT path with a stored profile: floor + tonal peaks are filled
    drainAll();
    proc.viz.feed(inGroup, std::span<const float>(input.data(), input.size()), std::span<const float>(output.data(), output.size()));
    drainAll();
    expect(frame.hasNoiseProfile, "profile frame reports the stored profile");
    expect(frame.hasNoiseFloorEstimate,
           "profile frame carries the floor estimate");
    expect(!frame.tonalPeaksHz.empty(),
           "tonal material must produce reported tonal peaks");
    expect(frame.tonalPeaksHz.size() <=
               NoiseRepellentAudioProcessor::kMaxTonalPeaks,
           "tonal peak list stays within its capacity");

    // Unlinked silent frame reports tonal peaks (visibility rule path)
    setParam(proc, "link_reduction", 0.0f);
    drainAll();
    for (int i = 0; i < 4; ++i) {
      proc.viz.pushSilentFrame(inGroup, proc.profiles);
    }
    drainAll();
    expect(frame.hasNoiseProfile, "frozen frame keeps the learned profile");
    expect(frame.hasNoiseFloorEstimate,
           "frozen frame carries the profile floor estimate");
    expect(!frame.isLinked, "unlinked flag must reach the frame");

    // Offset-unlinked silent frame (broadband linked, tonal offset free)
    setParam(proc, "link_reduction", 1.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    drainAll();
    for (int i = 0; i < 4; ++i) {
      proc.viz.pushSilentFrame(inGroup, proc.profiles);
    }
    drainAll();
    expect(frame.isLinked && !frame.isOffsetLinked,
           "offset link flag must reach the frame");

    // Linked silent frame skips the tonal peak query
    setParam(proc, "link_reduction", 1.0f);
    setParam(proc, "link_threshold_offset", 1.0f);
    drainAll();
    for (int i = 0; i < 4; ++i) {
      proc.viz.pushSilentFrame(inGroup, proc.profiles);
    }
    drainAll();
    expect(frame.isLinked && frame.isOffsetLinked,
           "link flags must reach the frame");

    // Frozen-profile path under lock contention: previous data retained
    {
      juce::SpinLock::ScopedLockType guard(proc.profiles.lock());
      for (int i = 0; i < 4; ++i) {
        proc.viz.pushSilentFrame(inGroup, proc.profiles);
      }
    }
    drainAll();
    proc.releaseResources();

    // Adaptive floor estimate without any stored profile, then FIFO overflow
    NoiseRepellentAudioProcessor adaptive;
    adaptive.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    VisualizationFeed::FeedInput inAdaptive{
        .group = adaptive.engine.group(),
        .profileLock = adaptive.profiles.lock(),
        .apvts = adaptive.getAPVTS(),
        .latency = 0,
        .isLearning = false,
        .adaptiveNoise = true,
        .transientIntensity = 0.0f,
        .transientProtectionActive = false};
    adaptive.viz.feed(inAdaptive, std::span<const float>(input.data(), input.size()), std::span<const float>(output.data(), output.size()));
    expect(adaptive.viz.getNextFrame(frame),
           "adaptive feed must push a frame");
    expect(!frame.hasNoiseProfile,
           "adaptive-only state has no stored profile");
    expect(frame.hasNoiseFloorEstimate,
           "adaptive noise still provides a floor estimate");

    // Fill the 16-slot ring buffer: full-frame pushes and silent pushes must
    // both tolerate a full FIFO
    for (int i = 0; i < 20; ++i) {
      adaptive.viz.feed(inAdaptive, std::span<const float>(input.data(), 1024), std::span<const float>(output.data(), 1024));
    }
    for (int i = 0; i < 72; ++i) {
      adaptive.viz.pushSilentFrame(inAdaptive, adaptive.profiles);
    }
    int drained = 0;
    while (adaptive.viz.getNextFrame(frame)) {
      ++drained;
    }
    expect(drained > 0, "overflowed ring buffer still yields frames");
    adaptive.releaseResources();

    // Feed with no prepared delay line at all
    NoiseRepellentAudioProcessor cold;
    VisualizationFeed::FeedInput inCold{
        .group = nullptr,
        .profileLock = cold.profiles.lock(),
        .apvts = cold.getAPVTS(),
        .latency = 0,
        .isLearning = false,
        .adaptiveNoise = false,
        .transientIntensity = 0.0f,
        .transientProtectionActive = false};
    cold.viz.feed(inCold, std::span<const float>(input.data(), input.size()), std::span<const float>());
    expect(cold.viz.getNextFrame(frame),
           "feed without a delay line must still push");
  }

  void testVizThroughProcessBlock() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    ScopedEditor editor(proc);
    pumpMessageLoop(20);
    expect(editor.editor != nullptr, "editor must open");

    juce::MidiBuffer midi;
    juce::AudioBuffer<float> buffer(2, 512);

    // Noise blocks accumulate into full FFT windows (past the latency
    // compensation delay so the aligned input carries real signal)
    for (int b = 0; b < 16; ++b) {
      generateNoiseBuffer(buffer, 0.05f, 100 + b);
      proc.processBlock(buffer, midi);
    }
    NoiseRepellentAudioProcessor::SpectralFrame frame;
    bool gotFrame = false;
    float loudestInputBin = -200.0f;
    while (proc.viz.getNextFrame(frame)) {
      gotFrame = true;
      loudestInputBin =
          std::max(loudestInputBin, frame.inputMagnitudeDB[100]);
    }
    expect(gotFrame, "editor-open processing must feed the visualization");
    expect(loudestInputBin > -120.0f,
           "noise input must show energy in the spectrum");

    // Silent blocks push the cheap silent frame path (the block must be
    // re-zeroed every iteration: processBlock overwrites it with output)
    for (int b = 0; b < 8; ++b) {
      buffer.clear();
      proc.processBlock(buffer, midi);
    }
    gotFrame = false;
    while (proc.viz.getNextFrame(frame)) {
      gotFrame = true;
    }
    expect(gotFrame, "silent blocks must push silent frames");
    expectEquals(frame.inputMagnitudeDB[0], -120.0f,
                 "silent frames floor the input spectrum");

    // Frame-size switch while the editor is open repaints through the rebuild
    setParam(proc, "frame_size", 3.0f);
    pumpMessageLoop(20);
    expectEquals(proc.getFrameSizeMs(), 64.0f,
                 "editor-open frame switch must rebuild");
    proc.releaseResources();
  }

  void testEngineParamCache() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    juce::AudioBuffer<float> buffer(2, 512);
    generateNoiseBuffer(buffer, 0.05f, 111);

    // Curve enabled with a mismatched size: the library rejects the load and
    // the engine keeps its previous parameters
    auto ep = proc.engine.buildParams(proc.getAPVTS());
    ep.p.reduction_curve_enabled = true;
    ep.p.reduction_curve_bias = nullptr;
    ep.p.reduction_curve_size = 0;
    proc.engine.run(buffer, ep);
    expectAllFinite(*this, buffer,
                    "rejected parameter load must not break processing");

    const auto base = proc.engine.buildParams(proc.getAPVTS());
    auto run = [&](EngineRuntime::EngineParams p, int seed) {
      generateNoiseBuffer(buffer, 0.05f, seed);
      proc.engine.run(buffer, p);
      expectAllFinite(*this, buffer, "engine run must stay finite");
    };
    int seed = 120;
    auto miss = [&](auto mutate) {
      auto changed = base;
      mutate(changed);
      run(changed, seed++);
      run(base, seed++); // restore the cache with the base snapshot
    };

    run(base, 112);
    run(base, 113); // steady state: cache fast path

    // Every compared field must be able to miss the cache on its own
    miss([](EngineRuntime::EngineParams& p) {
      p.p.learn_noise = SPECBLEACH_LEARN_ALL;
    });
    miss([](EngineRuntime::EngineParams& p) { p.p.residual_listen = true; });
    miss([](EngineRuntime::EngineParams& p) { p.p.reduction_gain *= 1.01f; });
    miss([](EngineRuntime::EngineParams& p) { p.p.smoothing_factor += 0.01f; });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.smoothing_mode = SPECBLEACH_SMOOTHING_NLM_2D;
    });
    miss([](EngineRuntime::EngineParams& p) { p.p.dftt_strength += 0.5f; });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.whitening_factor += 0.01f;
    });
    miss([](EngineRuntime::EngineParams& p) { p.p.adaptive_noise = true; });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.noise_estimation_method = SPECBLEACH_NOISE_ESTIMATION_SPP_MMSE;
    });
    miss([](EngineRuntime::EngineParams& p) { p.p.masking_depth += 0.01f; });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.suppression_strength = 0.5f;
    });
    miss([](EngineRuntime::EngineParams& p) { p.p.aggressiveness += 0.1f; });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.tonal_reduction_gain *= 1.01f;
    });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.transient_protection_enable = true;
    });
    miss([](EngineRuntime::EngineParams& p) { p.p.noise_profile_scale *= 1.01f; });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.tonal_noise_profile_scale *= 1.01f;
    });
    static float dummyBias = 0.0f;
    miss([](EngineRuntime::EngineParams& p) {
      p.p.reduction_curve_bias = &dummyBias;
    });
    miss([](EngineRuntime::EngineParams& p) {
      p.p.reduction_curve_size = 7; // curve disabled: size still compares
    });

    // Curve-enabled snapshots with missing bias / zero size exercise the
    // guard's short circuits before the library rejects the load
    auto noBias = base;
    noBias.curveEnabled = true;
    noBias.p.reduction_curve_enabled = true;
    noBias.p.reduction_curve_bias = nullptr;
    noBias.p.reduction_curve_size = 0;
    run(noBias, seed++);
    auto noSize = base;
    noSize.curveEnabled = true;
    noSize.p.reduction_curve_enabled = true;
    noSize.p.reduction_curve_bias = proc.curve.biasData();
    noSize.p.reduction_curve_size = 0;
    run(noSize, seed++);

    // Curve-enabled snapshots: load, cache, and curve-content mismatch
    proc.curve.setNodes({{0.3f, 2.0f}});
    proc.curve.refreshIfDirty(true, proc.engine.group(), 48000.0);
    auto withCurve = base;
    withCurve.curveEnabled = true;
    withCurve.p.reduction_curve_enabled = true;
    withCurve.p.reduction_curve_bias = proc.curve.biasData();
    withCurve.p.reduction_curve_size = proc.curve.biasSize();
    run(withCurve, seed++);
    run(withCurve, seed++); // curve content unchanged: cache hit
    proc.curve.setNodes({{0.3f, 2.0f}, {0.7f, -4.0f}});
    proc.curve.refreshIfDirty(true, proc.engine.group(), 48000.0);
    withCurve.p.reduction_curve_bias = proc.curve.biasData();
    withCurve.p.reduction_curve_size = proc.curve.biasSize();
    run(withCurve, seed++); // curve bytes changed: cache misses
    run(base, seed++);     // back to no curve: cache cleared
    proc.releaseResources();
  }

  void testFactoryAndEditor() {
    juce::AudioProcessor::setTypeOfNextNewPlugin(
        juce::AudioProcessor::wrapperType_Undefined);
    juce::AudioProcessor* raw = createPluginFilter();
    expect(raw != nullptr, "factory must create the processor");
    delete raw;

    NoiseRepellentAudioProcessor proc;
    {
      ScopedEditor editor(proc);
      pumpMessageLoop(20);
      expect(editor.editor != nullptr, "editor must be creatable");
    }
    expect(proc.getActiveEditor() == nullptr,
           "editor must detach cleanly after deletion");
    proc.releaseResources();
  }
};

ProcessorCoverageTest processorCoverageTest;

int main(int argc, char* argv[]) {
  juce::ignoreUnused(argc, argv);
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
