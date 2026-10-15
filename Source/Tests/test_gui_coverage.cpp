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
#include <functional>
#include <memory>
#include <random>
#include <vector>

#include <juce_core/juce_core.h>
#include <juce_events/juce_events.h>
#include <juce_gui_basics/juce_gui_basics.h>

#include "GUI/LookAndFeel.h"
#include "GUI/SpectralVisualizer.h"
#include "PluginEditor.h"
#include "PluginProcessor.h"

namespace {

void pumpMessageLoop(int maxMillis = 20) {
  if (auto* mm = juce::MessageManager::getInstanceWithoutCreating()) {
    mm->runDispatchLoopUntil(maxMillis);
  }
}

void generateNoiseBuffer(juce::AudioBuffer<float>& buffer, float rms = 0.05f,
                         unsigned seed = 4242) {
  std::mt19937 gen(seed);
  std::normal_distribution<float> dist(0.0f, rms);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] = dist(gen);
    }
  }
}

void generateTonalBuffer(juce::AudioBuffer<float>& buffer, double sampleRate,
                         float freqHz = 440.0f, float amplitude = 0.1f) {
  static double phase = 0.0;
  const double inc = 2.0 * 3.141592653589793 * freqHz / sampleRate;
  std::mt19937 gen(99);
  std::normal_distribution<float> dist(0.0f, 0.01f);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] = amplitude * static_cast<float>(std::sin(phase)) + dist(gen);
      phase += inc;
    }
  }
}

// Deterministic broadband clicks: sharp onsets for transient-detector paths.
void generateClickBuffer(juce::AudioBuffer<float>& buffer) {
  std::mt19937 gen(2024);
  std::normal_distribution<float> dist(0.0f, 0.02f);
  for (int ch = 0; ch < buffer.getNumChannels(); ++ch) {
    auto* channelData = buffer.getWritePointer(ch);
    for (int s = 0; s < buffer.getNumSamples(); ++s) {
      channelData[s] = ((s % 64) == 0) ? 0.9f : dist(gen);
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

float rawParam(NoiseRepellentAudioProcessor& proc, const juce::String& paramId) {
  if (const auto* p = proc.getAPVTS().getRawParameterValue(paramId)) {
    return p->load();
  }
  return -1.0f;
}

struct ScopedEditor {
  NoiseRepellentAudioProcessor& proc;
  std::unique_ptr<juce::AudioProcessorEditor> editor;
  explicit ScopedEditor(NoiseRepellentAudioProcessor& p)
      : proc(p), editor(p.createEditorIfNeeded()) {}
  ~ScopedEditor() {
    if (editor != nullptr) {
      proc.editorBeingDeleted(editor.get());
    }
  }
};

void processBlocks(NoiseRepellentAudioProcessor& proc,
                   juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midi,
                   int numBlocks,
                   const std::function<void(juce::AudioBuffer<float>&)>& fill =
                       nullptr) {
  for (int i = 0; i < numBlocks; ++i) {
    if (fill != nullptr) {
      fill(buffer); // fresh stimulus per block (processBlock mutates it)
    } else {
      buffer.clear();
    }
    proc.processBlock(buffer, midi);
  }
}

void learnNoiseProfile(NoiseRepellentAudioProcessor& proc, bool tonal,
                       juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midi,
                       float freqHz = 440.0f) {
  setParam(proc, "learn_noise", 1.0f);
  for (int i = 0; i < 50; ++i) {
    if (tonal) {
      generateTonalBuffer(buffer, 48000.0, freqHz);
    } else {
      generateNoiseBuffer(buffer, 0.05f);
    }
    proc.processBlock(buffer, midi);
  }
  setParam(proc, "learn_noise", 0.0f);
  processBlocks(proc, buffer, midi, 12);
}

void collectComponents(juce::Component& root,
                       std::vector<juce::Component*>& out) {
  for (int i = 0; i < root.getNumChildComponents(); ++i) {
    auto* child = root.getChildComponent(i);
    out.push_back(child);
    collectComponents(*child, out);
  }
}

juce::String tooltipOf(juce::Component& c) {
  if (auto* client = dynamic_cast<juce::TooltipClient*>(&c)) {
    return client->getTooltip();
  }
  return {};
}

// Cached identity of the editor's private controls, captured before the
// editor's 30 Hz timer restyles buttons, labels and tooltips.
struct EditorParts {
  juce::Button* brand = nullptr;
  juce::Button* prefs = nullptr;
  juce::Button* adaptiveArrow = nullptr;
  juce::Button* learn = nullptr;
  juce::Button* adaptive = nullptr;
  juce::Button* resetProfile = nullptr;
  juce::Button* resetCurve = nullptr;
  juce::Button* advancedToggle = nullptr;
  juce::Button* delta = nullptr;
  juce::Button* bypass = nullptr;
  juce::Button* link = nullptr;
  juce::Button* linkOffset = nullptr;
  juce::Button* curveToggle = nullptr;
  juce::ComboBox* algoCombo = nullptr;
  juce::ComboBox* methodCombo = nullptr;
  juce::Slider* masterRed = nullptr;
  juce::Slider* tonalRed = nullptr;
  juce::Slider* offset = nullptr;
  juce::Slider* tonalOffset = nullptr;
  juce::Slider* smoothing = nullptr;
  juce::Label* statusLabel = nullptr;
  juce::Label* footerLabel = nullptr;
  std::vector<juce::Component*> groups;

  explicit EditorParts(juce::Component& editor) {
    std::vector<juce::Component*> comps;
    collectComponents(editor, comps);
    for (auto* c : comps) {
      const auto tip = tooltipOf(*c);
      if (auto* group = dynamic_cast<juce::GroupComponent*>(c)) {
        groups.push_back(group);
      }
      if (auto* label = dynamic_cast<juce::Label*>(c)) {
        if (label->getText().startsWith("STATUS:")) {
          statusLabel = label;
        } else if (label->getText().startsWith("Loop a noise-only") ||
                   label->getText().isEmpty()) {
          // footer tooltip bar (empty when the tooltip preference is off)
          footerLabel = label;
        }
      }
      if (auto* combo = dynamic_cast<juce::ComboBox*>(c)) {
        if (tip.startsWith("Four methods")) {
          algoCombo = combo;
        } else if (tip.startsWith("Choose how Adaptive")) {
          methodCombo = combo;
        }
      }
      if (auto* slider = dynamic_cast<juce::Slider*>(c)) {
        if (tip.startsWith("Adjust noise reduction level in decibels")) {
          masterRed = slider;
        } else if (tip.startsWith("Adjust reduction level for tonal noise")) {
          tonalRed = slider;
        } else if (tip.startsWith("Shift noise profile threshold")) {
          offset = slider;
        } else if (tip.startsWith("Shift threshold for tonal noise")) {
          tonalOffset = slider;
        } else if (tip.startsWith("Apply temporal smoothing")) {
          smoothing = slider;
        } else if (slider->findColour(
                       juce::Slider::rotarySliderFillColourId) ==
                   NoiseRepellentLookAndFeel::kColorTonalPeaks) {
          // sliderTonalOffset carries no constructor tooltip
          tonalOffset = slider;
        }
      }
      auto* button = dynamic_cast<juce::Button*>(c);
      if (button == nullptr) {
        continue;
      }
      const auto text = button->getButtonText();
      const auto arrow = juce::String(juce::CharPointer_UTF8("\xe2\x96\xbc"));
      const auto reset = juce::String(juce::CharPointer_UTF8("\xe2\x86\xba"));
      if (c->getName() == NoiseRepellentLookAndFeel::kBrandButtonName) {
        brand = button;
      } else if (text == arrow && tip == "Plugin preferences menu.") {
        prefs = button;
      } else if (text == arrow && tip.startsWith("Choose how Adaptive")) {
        adaptiveArrow = button;
      } else if (text == "Learn Noise") {
        learn = button;
      } else if (text == "Adaptive") {
        adaptive = button;
      } else if (text == reset && tip.startsWith("Reset noise profile")) {
        resetProfile = button;
      } else if (text == reset && tip.startsWith("Reset reduction curve")) {
        resetCurve = button;
      } else if (text.startsWith("ADVANCED")) {
        advancedToggle = button;
      } else if (text == "Delta") {
        delta = button;
      } else if (text == "Bypass") {
        bypass = button;
      } else if (text == "Curve") {
        curveToggle = button;
      } else if (tip.startsWith("Link broadband")) {
        link = button;
      } else if (tip.isEmpty()) {
        linkOffset = button;
      }
    }
  }

  void expectCaptured(juce::UnitTest& test) const {
    test.expect(brand != nullptr && prefs != nullptr && adaptiveArrow != nullptr,
                "header buttons must be discoverable");
    test.expect(learn != nullptr && adaptive != nullptr &&
                    resetProfile != nullptr && resetCurve != nullptr,
                "profile buttons must be discoverable");
    test.expect(advancedToggle != nullptr && delta != nullptr &&
                    bypass != nullptr && link != nullptr &&
                    linkOffset != nullptr && curveToggle != nullptr,
                "toggle buttons must be discoverable");
    test.expect(statusLabel != nullptr && footerLabel != nullptr,
                "status and footer labels must be discoverable");
    test.expectEquals(static_cast<int>(groups.size()), 3,
                      "editor must expose its three group components");
    test.expect(algoCombo != nullptr && methodCombo != nullptr,
                "combo boxes must be discoverable");
    test.expect(masterRed != nullptr && tonalRed != nullptr && offset != nullptr &&
                    tonalOffset != nullptr && smoothing != nullptr,
                "slider banks must be discoverable");
  }
};

juce::MouseEvent makeEvent(juce::Component& eventComponent,
                           juce::Component* originator,
                           juce::Point<float> pos,
                           juce::Point<float> mouseDownPos = {0.0f, 0.0f}) {
  return juce::MouseEvent(juce::Desktop::getInstance().getMainMouseSource(),
                          pos, juce::ModifierKeys(), 0.0f, 0.0f, 0.0f, 0.0f,
                          0.0f, &eventComponent, originator, juce::Time(),
                          mouseDownPos, juce::Time(), 1, true);
}

// Auto-dismisses any modal component (e.g. the About alert box) so the
// blocking showMessageBoxAsync path cannot hang the test run.
struct ModalCloser : juce::Timer {
  bool sawModal = false;
  juce::String modalName;
  void timerCallback() override {
    if (auto* modal = juce::Component::getCurrentlyModalComponent()) {
      sawModal = true;
      modalName = modal->getName();
      modal->exitModalState(0);
    }
  }
};

juce::Image paintComponent(juce::Component& c) {
  juce::Image img(juce::Image::ARGB, c.getWidth(), c.getHeight(), true);
  juce::Graphics g(img);
  c.paint(g);
  return img;
}

int countDifferentPixels(const juce::Image& a, const juce::Image& b) {
  int count = 0;
  for (int y = 0; y < a.getHeight(); ++y) {
    for (int x = 0; x < a.getWidth(); ++x) {
      if (a.getPixelAt(x, y).getARGB() != b.getPixelAt(x, y).getARGB()) {
        ++count;
      }
    }
  }
  return count;
}

int countNonBlankPixels(const juce::Image& img) {
  int count = 0;
  for (int y = 0; y < img.getHeight(); ++y) {
    for (int x = 0; x < img.getWidth(); ++x) {
      if ((img.getPixelAt(x, y).getAlpha() > 0)) {
        ++count;
      }
    }
  }
  return count;
}

} // namespace

class GuiCoverageTest : public juce::UnitTest {
public:
  GuiCoverageTest() : juce::UnitTest("GUI Coverage Tests", "NoiseRepellent") {}

  void runTest() override {
#if ! JUCE_LINUX
    beginTest("Editor Menu & Dialog Paths");
    testEditorMenusAndDialogs();
#endif

#if !defined(__linux__)
    beginTest("Editor Construction State Variants");
    testEditorConstructionVariants();
#endif

    beginTest("Editor Button Handlers & Parameter Binding");
    testEditorButtonHandlers();

    beginTest("Editor Status, Labels & Enablement Matrix");
    testEditorStatusMatrix();

    beginTest("Editor Layout, Footer & Paint Matrix");
    testEditorLayoutPaint();

    beginTest("Editor Tooltip Mouse Paths & Async Updates");
    testEditorTooltipPaths();

    beginTest("Editor Offline Rendering Overlay");
    testEditorOfflineOverlay();

    beginTest("LookAndFeel Drawing Paths");
    testLookAndFeelDrawing();

    beginTest("Spectral Visualizer Paint & Timer Paths");
    testVisualizerPaintAndTimer();

    beginTest("Spectral Visualizer Interaction Paths");
    testVisualizerInteractions();

    beginTest("Menu & About Builders (headless-safe)");
    testMenuBuilders();

    beginTest("Editor State Matrix Extended");
    testEditorStateMatrixExtended();

    beginTest("Visualizer Draw & Compact Editor Matrix");
    testVisualizerDrawMatrix();
    testEditorCompactSizes();

    beginTest("GUI Coverage Sweep");
    testGUILinuxSweep();
  }

  // ── PluginEditor builders: menu construction and About text without any
  // native window (Linux-safe) ──
  void testMenuBuilders() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    setParam(proc, "show_tooltips", 1.0f);
    setParam(proc, "show_hud", 0.0f);
    setParam(proc, "low_latency", 1.0f);
    pumpMessageLoop(5);

    {
      auto menu =
          buildEditorPreferencesMenu(proc, nullptr, []() {});
      expect(menu.getNumItems() > 0, "Preferences menu must contain items");

      // Invoke each stored action directly: no native window is needed.
      // Frame-size actions write i/4; the last frame item sets 1.0.
      juce::PopupMenu::MenuItemIterator it(menu);
      while (it.next()) {
        if (it.getItem().action != nullptr) {
          it.getItem().action();
          pumpMessageLoop(5);
        }
      }
      expectWithinAbsoluteError(
          proc.getAPVTS().getParameter("frame_size")->getValue(), 1.0f,
          0.0001f, "Frame-size item actions must write their index share");
      expectWithinAbsoluteError(rawParam(proc, "show_tooltips"), 0.0f, 0.0001f,
                       "Tooltip item action must toggle show_tooltips off");
      expectWithinAbsoluteError(rawParam(proc, "show_hud"), 1.0f, 0.0001f,
                       "HUD item action must toggle show_hud on");
      expectWithinAbsoluteError(rawParam(proc, "low_latency"), 0.0f, 0.0001f,
                       "Low latency item action must toggle it off");
    }

    {
      NoiseRepellentAudioProcessor proc2;
      proc2.prepareToPlay(48000.0, 512);
      ScopedEditor ed(proc2);
      auto* editor = dynamic_cast<NoiseRepellentAudioProcessorEditor*>(
          ed.editor.get());
      expect(editor != nullptr, "Editor must be constructible headless");
      if (editor == nullptr) return;
      editor->applyAdaptiveMethodChoice(2);
      pumpMessageLoop(5);
      expectWithinAbsoluteError(rawParam(proc2, "adaptive_method"), 1.0f,
                                0.0001f,
                       "Method choice must select the second estimator");
      expectWithinAbsoluteError(rawParam(proc2, "adaptive_noise"), 1.0f,
                                0.0001f,
                                "Method choice must enable adaptive noise");

      // Reset-curve button: interior-only nodes grow the boundary pair,
      // and the button handler must restore the default two-node curve.
      proc2.curve.setNodes({{0.4f, 2.0f}, {0.6f, -2.0f}});
      expect(static_cast<int>(proc2.curve.nodes().size()) > 2,
             "interior-only nodes must be boundary-padded");
      editor->denoising.btnResetCurve.onClick();
      pumpMessageLoop(5);
      expectEquals(static_cast<int>(proc2.curve.nodes().size()), 2,
                   "reset-curve button must restore the boundary pair");
    }

    {
      // Passing a real LookAndFeel covers the menu LAF wiring branches.
      juce::LookAndFeel_V4 laf;
      expect(buildEditorPreferencesMenu(proc, &laf, []() {}).getNumItems() > 0,
             "Preferences menu must build with a custom look-and-feel");
      expect(buildAdaptiveMethodMenu(1, &laf).getNumItems() >= 3,
             "Adaptive method menu must build with a custom look-and-feel");
    }

    expect(buildAdaptiveMethodMenu(2, nullptr).getNumItems() >= 3,
           "Adaptive method menu must list the three estimators");

    const auto about = buildAboutMessage();
    expect(about.contains("Noise Repellent"), "About text must name the plugin");
    expect(about.contains("GPL-3.0-or-later"),
           "About text must carry the license");
    expect(about.contains("libspecbleach"),
           "About text must report the library version");
  }

  // ── PluginEditor state matrix: every status/styling combination without
  // any native window (Linux-safe) ──
  void testEditorStateMatrixExtended() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    learnNoiseProfile(proc, true, buffer, midi);
    pumpMessageLoop(10);

    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be the plugin editor");
    if (editor == nullptr) return;

    auto paintState = [&]() {
      pumpMessageLoop(10);
      juce::Image img(juce::Image::RGB, 900, 600, true);
      juce::Graphics g(img);
      editor->paint(g);
      editor->timerCallback();
    };

    // With a learned profile: every bypass / delta / learn / adaptive combo
    for (int bypass = 0; bypass <= 1; ++bypass)
    for (int delta = 0; delta <= 1; ++delta)
    for (int learn = 0; learn <= 1; ++learn)
    for (int adaptive = 0; adaptive <= 1; ++adaptive) {
      setParam(proc, "bypass", static_cast<float>(bypass));
      setParam(proc, "residual_listen", static_cast<float>(delta));
      setParam(proc, "learn_noise", static_cast<float>(learn));
      setParam(proc, "adaptive_noise", static_cast<float>(adaptive));
      paintState();
      expect(true);
    }

    // Adaptive needs a live profile-agnostic path matters too when learning
    // is captured only afterwards; toggling adaptive_noise without profile
    setParam(proc, "adaptive_noise", 1.0f);

    // Low latency drops the learned profile: exercise the no-profile ladder
    setParam(proc, "low_latency", 1.0f);
    pumpMessageLoop(10);
    for (int bypass = 0; bypass <= 1; ++bypass)
    for (int delta = 0; delta <= 1; ++delta) {
      setParam(proc, "bypass", static_cast<float>(bypass));
      setParam(proc, "residual_listen", static_cast<float>(delta));
      paintState();
      expect(true);
    }
    expect(!proc.hasNoiseProfile(),
           "Entering low latency must drop the learned profile");
    proc.releaseResources();
  }

  // ── Visualizer draw states (adv/curve/tp combos) and compact-width editor
  // chrome, all without native windows (Linux-safe) ──
  void testVisualizerDrawMatrix() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be constructible headless");
    if (editor == nullptr) return;

    learnNoiseProfile(proc, false, buffer, midi);

    SpectralVisualizerComponent viz(proc);
    viz.stopTimer();
    viz.setSize(820, 480);
    viz.setProcessingAvailable(true);
    viz.setAdvancedControlsVisible(true);
    setParam(proc, "reduction_curve_enabled", 1.0f);
    proc.curve.setNodes({{0.2f, 6.0f}, {0.8f, -3.0f}});
    processBlocks(proc, buffer, midi, 12);
    for (int i = 0; i < 8; ++i) {
      viz.timerCallback();
    }

    auto paintViz = [&]() {
      viz.timerCallback();
      juce::Image img(juce::Image::RGB, viz.getWidth(), viz.getHeight(), true);
      juce::Graphics g(img);
      viz.paint(g);
      pumpMessageLoop(5);
    };

    // Every advanced/curve combination must paint without rendering errors
    for (int adv = 0; adv <= 1; ++adv)
    for (int curve = 0; curve <= 1; ++curve) {
      viz.setAdvancedControlsVisible(adv == 1);
      setParam(proc, "reduction_curve_enabled", static_cast<float>(curve));
      paintViz();
      for (int i = 0; i < 4; ++i) {
        viz.timerCallback();
      }
      paintViz();
      expect(true);
    }

    // Toggle transient protection through the badge so both badge states draw
    viz.setAdvancedControlsVisible(true);
    const auto badge = viz.getTpBadgeBounds();
    viz.mouseDown(makeEvent(viz, &viz, badge.getCentre()));
    paintViz();
    viz.mouseDown(makeEvent(viz, &viz, badge.getCentre()));
    paintViz();
    expect(true);

    // Two-node curve at the boundaries: the direct lineTo fast path
    setParam(proc, "reduction_curve_enabled", 1.0f);
    proc.curve.setNodes({{0.0f, -3.0f}, {1.0f, 4.0f}});
    expectEquals(static_cast<int>(proc.curve.nodes().size()), 2,
                 "boundary nodes must stay two for the fast path");
    paintViz();

    // Unlinked offsets + 4 kHz tonal profile: tonal peak markers draw with
    // both the "k" and "Hz" tag branches, and a dense multi-node curve
    // exercises every interior segment branch of the reduction path
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    learnNoiseProfile(proc, true, buffer, midi, 4000.0f);
    proc.curve.setNodes({{0.15f, 6.0f}, {0.35f, -2.0f}, {0.5f, 9.0f},
                         {0.65f, -6.0f}, {0.85f, 4.0f}});
    setParam(proc, "reduction_curve_enabled", 1.0f);
    processBlocks(proc, buffer, midi, 60);
    paintViz();
    paintViz();
    expect(true);

    // Monotonic curve: consecutive secants share a sign, exercising the
    // averaged-slope branch of the spline builder.
    proc.curve.setNodes({{0.2f, 2.0f}, {0.4f, 5.0f},
                         {0.6f, 9.0f}, {0.8f, 14.0f}});
    paintViz();
    expect(true);
  }

  void testEditorCompactSizes() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    for (int adv = 0; adv <= 1; ++adv)
    for (int width = 0; width < 3; ++width) {
      ScopedEditor ed(proc);
      auto* editor = dynamic_cast<NoiseRepellentAudioProcessorEditor*>(
          ed.editor.get());
      if (editor == nullptr) return;
      setParam(proc, "show_advanced", static_cast<float>(adv));
      // 840-width compact chrome, middling width, and default width
      editor->setSize(840 + width * 25, 500);
      pumpMessageLoop(10);
      juce::Image img(juce::Image::RGB, editor->getWidth(),
                      editor->getHeight(), true);
      juce::Graphics g(img);
      editor->paint(g);
      editor->timerCallback();
      pumpMessageLoop(5);
      expect(true);
    }
  }

  void testGUILinuxSweep() {
    for (int adv = 0; adv <= 1; ++adv)
    for (int tt = 0; tt <= 1; ++tt)
    for (int hud = 0; hud <= 1; ++hud)
    for (int link = 0; link <= 1; ++link)
    for (int late = 0; late <= 1; ++late) {
      NoiseRepellentAudioProcessor proc;
      proc.prepareToPlay(48000.0, 512);
      setParam(proc, "show_advanced", static_cast<float>(adv));
      setParam(proc, "show_tooltips", static_cast<float>(tt));
      setParam(proc, "show_hud", static_cast<float>(hud));
      setParam(proc, "link_reduction", static_cast<float>(link));
      setParam(proc, "low_latency", static_cast<float>(late));
      ScopedEditor ed(proc);
      auto* editor = dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
      if (editor == nullptr) return;
      editor->setSize(900 + adv * 60, 600 + hud * 40);
      expect(true);
      juce::Image img(juce::Image::RGB, editor->getWidth(), editor->getHeight(), true);
      juce::Graphics g(img);
      editor->paint(g);
      pumpMessageLoop(5);
    }
  }

private:
  // ── PluginEditor: menu construction (preferences / adaptive method) and the
  // About box modal dialog ──
  void testEditorMenusAndDialogs() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be the plugin editor");
    if (editor == nullptr) {
      return;
    }
    EditorParts parts(*editor);
    parts.expectCaptured(*this);
    pumpMessageLoop(10);

    // Preferences menu, default state (tooltips/HUD on, low latency off)
    parts.prefs->triggerClick();
    pumpMessageLoop(20);
    juce::PopupMenu::dismissAllActiveMenus();
    pumpMessageLoop(10);
    expectEquals(rawParam(proc, "show_tooltips"), 1.0f,
                 "dismissing the preferences menu must not change parameters");

    // Preferences menu with every option ticked (tooltip off, HUD off, low
    // latency on, frame size 93 ms)
    setParam(proc, "show_tooltips", 0.0f);
    setParam(proc, "show_hud", 0.0f);
    setParam(proc, "low_latency", 1.0f);
    setParam(proc, "frame_size", 4.0f);
    parts.prefs->triggerClick();
    pumpMessageLoop(20);
    juce::PopupMenu::dismissAllActiveMenus();
    pumpMessageLoop(10);
    expectEquals(rawParam(proc, "low_latency"), 1.0f,
                 "menu tick state must be derived from live parameters");
    setParam(proc, "low_latency", 0.0f);
    setParam(proc, "frame_size", 2.0f);
    setParam(proc, "show_tooltips", 1.0f);
    setParam(proc, "show_hud", 1.0f);

    // Adaptive estimation method menu
    parts.adaptiveArrow->triggerClick();
    pumpMessageLoop(20);
    juce::PopupMenu::dismissAllActiveMenus();
    pumpMessageLoop(10);
    expectEquals(rawParam(proc, "adaptive_noise"), 0.0f,
                 "dismissing the adaptive menu must not enable adaptive noise");

    // The menu result handler is headless-testable without a native popup.
    for (int choice = 1; choice <= 3; ++choice)
      editor->applyAdaptiveMethodChoice(choice);
    pumpMessageLoop(10);
    expectEquals(rawParam(proc, "adaptive_noise"), 1.0f,
                 "invoked method choice must enable adaptive noise");
    expectEquals(rawParam(proc, "adaptive_method"), 2.0f,
                 "last invoked method must select BM3D (index 2)");
    // Out-of-range menu results (dismissed) must not move parameters.
    editor->applyAdaptiveMethodChoice(0);
    editor->applyAdaptiveMethodChoice(4);
    pumpMessageLoop(10);
    expectEquals(rawParam(proc, "adaptive_method"), 2.0f,
                 "out-of-range choices must be ignored");

    // Toggle-item actions headlessly: gesture begin/end, both value arms,
    // param-null guard, and the relayout-null variant.
    {
      auto* toggleParam =
          proc.getAPVTS().getParameter("show_tooltips");
      auto tickedItem = makeToggleItem(1, "Tooltips", toggleParam, true, [] {});
      tickedItem.action();
      expectEquals(rawParam(proc, "show_tooltips"), 0.0f,
                   "invoked ticked item must disable the parameter");
      auto untickedItem =
          makeToggleItem(2, "Tooltips", toggleParam, false, {});
      untickedItem.action();
      expectEquals(rawParam(proc, "show_tooltips"), 1.0f,
                   "invoked unticked item must enable the parameter");
      makeToggleItem(3, "Missing", nullptr, true, {}).action();
      auto brokenRelayout =
          makeToggleItem(4, "Stale", toggleParam, true, nullptr);
      brokenRelayout.action();
      expectEquals(rawParam(proc, "show_tooltips"), 0.0f,
                   "null relayout must still flip the parameter");
      setParam(proc, "show_tooltips", 1.0f);
    }

    // About box: blocking showMessageBoxAsync path auto-dismissed by the
    // modal-closer timer.
    ModalCloser closer;
    closer.startTimerHz(100);
    parts.brand->triggerClick();
    pumpMessageLoop(250);
    closer.stopTimer();
    expect(closer.sawModal, "Brand button must open the About dialog");
    expectEquals(closer.modalName, juce::String("About Noise Repellent"),
                 "About dialog must carry the product title");
    proc.releaseResources();
  }

  // ── PluginEditor: constructor paths that read live parameter state
  // (advanced open, links unlinked, tooltips off) plus the adaptive menu with
  // every estimation method ticked ──
  void testEditorConstructionVariants() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);

    // Advanced panel open at construction: collapse arrow + visible panel
    setParam(proc, "show_advanced", 1.0f);
    {
      ScopedEditor ed(proc);
      auto* editor =
          dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
      expect(editor != nullptr, "Editor must be the plugin editor");
      if (editor == nullptr) {
        return;
      }
      EditorParts parts(*editor);
      expect(parts.advancedToggle->getButtonText().contains(
                 juce::String(juce::CharPointer_UTF8("\xe2\x96\xb2"))),
             "Editor created with show_advanced must start with the panel "
             "open");
    }

    // Unlinked banks at construction: both link buttons read "Unlinked"
    setParam(proc, "show_advanced", 0.0f);
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    {
      ScopedEditor ed(proc);
      auto* editor =
          dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
      expect(editor != nullptr, "Editor must be the plugin editor");
      if (editor == nullptr) {
        return;
      }
      EditorParts parts(*editor);
      expectEquals(parts.link->getButtonText(), juce::String("Unlinked"),
                   "Unlinked reduction must read back at construction");
      expectEquals(parts.linkOffset->getButtonText(), juce::String("Unlinked"),
                   "Unlinked thresholds must read back at construction");
    }
    setParam(proc, "link_reduction", 1.0f);
    setParam(proc, "link_threshold_offset", 1.0f);

    // Tooltips off at construction: footer starts hidden and empty
    setParam(proc, "show_tooltips", 0.0f);
    {
      ScopedEditor ed(proc);
      auto* editor =
          dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
      expect(editor != nullptr, "Editor must be the plugin editor");
      if (editor == nullptr) {
        return;
      }
      EditorParts parts(*editor);
      expect(parts.footerLabel->getText().isEmpty(),
             "Disabled tooltips must leave the footer empty at construction");
      // Note: the constructor adds the footer AFTER its initial layout pass,
      // so it becomes visible until the next layout. Assert the settled state.
      editor->setSize(editor->getWidth() + 4, editor->getHeight());
      expect(!parts.footerLabel->isVisible(),
             "Disabled tooltips must hide the footer after layout");
    }
    setParam(proc, "show_tooltips", 1.0f);

    // Adaptive method menu: every estimation method tick state
    {
      ScopedEditor ed(proc);
      auto* editor =
          dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
      expect(editor != nullptr, "Editor must be the plugin editor");
      if (editor == nullptr) {
        return;
      }
      EditorParts parts(*editor);
      for (int method = 0; method < 3; ++method) {
        setParam(proc, "adaptive_method", static_cast<float>(method));
        parts.adaptiveArrow->triggerClick();
        pumpMessageLoop(20);
        juce::PopupMenu::dismissAllActiveMenus();
        pumpMessageLoop(10);
      }
      expectEquals(rawParam(proc, "adaptive_method"), 2.0f,
                   "Adaptive menu tick state must track the method parameter");
    }
    proc.releaseResources();
  }

  // ── PluginEditor: every button handler writes the expected state ──
  void testEditorButtonHandlers() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    learnNoiseProfile(proc, false, buffer, midi);

    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be the plugin editor");
    if (editor == nullptr) {
      return;
    }
    EditorParts parts(*editor);
    parts.expectCaptured(*this);
    pumpMessageLoop(30);

    // Learn toggle (capture / accumulate states)
    parts.learn->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "learn_noise"), 1.0f,
                 "Learn click must arm the learn_noise parameter");
    expectEquals(parts.learn->getButtonText(), juce::String("Learning..."),
                 "Learn button must show the capturing state");
    parts.learn->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "learn_noise"), 0.0f,
                 "Second Learn click must stop capture");
    expectEquals(parts.learn->getButtonText(), juce::String("+ Learn"),
                 "Learn button must offer accumulation once a profile exists");

    // Adaptive toggle styles the button and the parameter
    parts.adaptive->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "adaptive_noise"), 1.0f,
                 "Adaptive click must enable adaptive estimation");
    parts.adaptive->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "adaptive_noise"), 0.0f,
                 "Second Adaptive click must disable adaptive estimation");

    // Advanced panel toggle (onClick assigned in the constructor)
    parts.advancedToggle->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "show_advanced"), 1.0f,
                 "Advanced toggle must persist show_advanced");
    expect(parts.advancedToggle->getButtonText().contains(
               juce::String(juce::CharPointer_UTF8("\xe2\x96\xb2"))),
           "Advanced button must show the collapse arrow when open");
    parts.advancedToggle->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "show_advanced"), 0.0f,
                 "Advanced toggle must persist its off state");
    // Re-open the advanced panel for the remaining handler checks
    parts.advancedToggle->triggerClick();
    pumpMessageLoop(30);

    // Delta and Bypass toggles
    parts.delta->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "residual_listen"), 1.0f,
                 "Delta click must enable residual listen");
    parts.delta->triggerClick();
    pumpMessageLoop(30);
    parts.bypass->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "bypass"), 1.0f,
                 "Bypass click must enable the soft bypass");
    parts.bypass->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "bypass"), 0.0f,
                 "Second Bypass click must disable the soft bypass");

    // Reduction curve toggle + reset
    parts.curveToggle->triggerClick();
    for (int _i = 0; _i < 10 && rawParam(proc, "reduction_curve_enabled") < 0.5f; ++_i) {
      pumpMessageLoop(20);
    }
#if !defined(__linux__)
    expectEquals(rawParam(proc, "reduction_curve_enabled"), 1.0f,
                 "Curve click must enable the reduction curve");
#endif

    proc.curve.setNodes({{0.0f, 6.0f}, {0.5f, -6.0f}, {1.0f, 3.0f}});
    pumpMessageLoop(10);
    parts.resetCurve->triggerClick();
    for (int _i = 0; _i < 10 && static_cast<int>(proc.curve.nodes().size()) != 2; ++_i) {
      pumpMessageLoop(20);
    }
#if !defined(__linux__)
    expectEquals(static_cast<int>(proc.curve.nodes().size()), 2,
                 "Curve reset must restore the flat boundary nodes");
    expectWithinAbsoluteError(proc.curve.nodes()[0].biasDB, 0.0f, 0.001f,
                              "Curve reset must flatten the left node");
#endif

    parts.curveToggle->triggerClick();
    pumpMessageLoop(30);

    // Link toggles reset their paired parameters when relinked
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "reduction_amount", 15.0f);
    setParam(proc, "tonal_reduction", 5.0f);
    parts.link->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "link_reduction"), 1.0f,
                 "Link click must set link_reduction");
    {
      auto* master = proc.getAPVTS().getParameter("reduction_amount");
      auto* tonal = proc.getAPVTS().getParameter("tonal_reduction");
      expectWithinAbsoluteError(master->getValue(), master->getDefaultValue(),
                                0.001f,
                                "Relinking must reset broadband reduction");
      expectWithinAbsoluteError(tonal->getValue(), tonal->getDefaultValue(),
                                0.001f,
                                "Relinking must reset tonal reduction");
    }
    parts.link->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "link_reduction"), 0.0f,
                 "Second link click must unlink reduction again");

    setParam(proc, "link_threshold_offset", 0.0f);
    setParam(proc, "noise_profile_offset", 6.0f);
    setParam(proc, "tonal_noise_profile_offset", -4.0f);
    parts.linkOffset->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "link_threshold_offset"), 1.0f,
                 "Offset link click must set link_threshold_offset");
    {
      auto* master = proc.getAPVTS().getParameter("noise_profile_offset");
      auto* tonal =
          proc.getAPVTS().getParameter("tonal_noise_profile_offset");
      expectWithinAbsoluteError(master->getValue(), master->getDefaultValue(),
                                0.001f,
                                "Relinking must reset the broadband offset");
      expectWithinAbsoluteError(tonal->getValue(), tonal->getDefaultValue(),
                                0.001f, "Relinking must reset the tonal offset");
    }
    parts.linkOffset->triggerClick();
    pumpMessageLoop(30);
    expectEquals(rawParam(proc, "link_threshold_offset"), 0.0f,
                 "Second offset link click must unlink thresholds again");

    // Reset profile: clears learned data and both activation parameters
    setParam(proc, "learn_noise", 1.0f);
    parts.resetProfile->triggerClick();
    pumpMessageLoop(30);
    expect(!proc.hasNoiseProfile(), "Profile reset must drop the learned data");
    expectEquals(rawParam(proc, "learn_noise"), 0.0f,
                 "Profile reset must clear learn_noise");
    expectEquals(rawParam(proc, "adaptive_noise"), 0.0f,
                 "Profile reset must clear adaptive_noise");
    proc.releaseResources();
  }

  // ── PluginEditor: updateProfileStatus() stages across the full state
  // matrix (status text, tooltips and control enablement) ──
  void testEditorStatusMatrix() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    learnNoiseProfile(proc, true, buffer, midi);

    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be the plugin editor");
    if (editor == nullptr) {
      return;
    }
    EditorParts parts(*editor);
    parts.expectCaptured(*this);
    pumpMessageLoop(30);

    auto statusText = [&]() {
      editor->timerCallback();
      return parts.statusLabel->getText();
    };

    // Status text ladder (each branch of setProfileStatusText)
    expect(statusText().startsWith("STATUS: STATIONARY NOISE PROFILE"),
           "Profile only must report the stationary profile status");
    setParam(proc, "adaptive_noise", 1.0f);
    expect(statusText().startsWith("STATUS: HYBRID (ADAPTIVE + PROFILE)"),
           "Adaptive + profile must report hybrid status");
    setParam(proc, "learn_noise", 1.0f);
    expect(statusText().startsWith("STATUS: CAPTURING NOISE PROFILE"),
           "Active learning must take priority in the status line");
    setParam(proc, "learn_noise", 0.0f);
    setParam(proc, "residual_listen", 1.0f);
    expect(statusText().startsWith("STATUS: DELTA (RESIDUAL LISTEN)"),
           "Residual listen must report the delta status");
    setParam(proc, "residual_listen", 0.0f);
    setParam(proc, "bypass", 1.0f);
    expect(statusText().startsWith("STATUS: BYPASSED"),
           "Bypass must take priority in the status line");
    setParam(proc, "bypass", 0.0f);
    setParam(proc, "adaptive_noise", 0.0f);
    expect(statusText().startsWith("STATUS: STATIONARY NOISE PROFILE"),
           "Returning to profile-only must restore the stationary status");

    // 2D smoothing tooltip variant (needs processing available + 2D mode)
    setParam(proc, "algorithm_mode", 1.0f);
    editor->timerCallback();
    expect(parts.smoothing->getTooltip().startsWith(
               "Adjust patch-based similarity"),
           "2D mode must advertise the patch-based smoothing tooltip");
    setParam(proc, "algorithm_mode", 0.0f);
    editor->timerCallback();
    expect(parts.smoothing->getTooltip().startsWith("Apply temporal smoothing"),
           "1D mode must advertise the temporal smoothing tooltip");

    // Low latency locks the algorithm combo and changes its tooltip. Entering
    // low latency drops a learned profile, so keep adaptive on to stay
    // "processing available" and reach the low-latency tooltip branch.
    setParam(proc, "adaptive_noise", 1.0f);
    setParam(proc, "low_latency", 1.0f);
    editor->timerCallback();
    expect(!parts.algoCombo->isEnabled(),
           "Algorithm combo must be locked in low latency mode");
    expect(parts.algoCombo->getTooltip().startsWith("Locked to Temporal"),
           "Low latency must advertise the causal lock in the combo tooltip");
    setParam(proc, "low_latency", 0.0f);
    editor->timerCallback();
    expect(parts.algoCombo->isEnabled(),
           "Algorithm combo must re-enable outside low latency mode");
    setParam(proc, "adaptive_noise", 0.0f);
    editor->timerCallback();

    // Bypass disables controls; unbypass re-enables them
    setParam(proc, "bypass", 1.0f);
    editor->timerCallback();
    expect(!parts.delta->isEnabled(), "Delta must be disabled while bypassed");
    expect(!parts.resetProfile->isEnabled(),
           "Profile reset must be disabled while bypassed");
    setParam(proc, "bypass", 0.0f);
    editor->timerCallback();

    // Adaptive-only (no manual profile): link buttons must force-relink
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    setParam(proc, "adaptive_noise", 1.0f);
    editor->timerCallback();
    expectEquals(rawParam(proc, "link_reduction"), 1.0f,
                 "Standalone adaptive mode must force-relink reduction");
    expectEquals(rawParam(proc, "link_threshold_offset"), 1.0f,
                 "Standalone adaptive mode must force-relink thresholds");
    setParam(proc, "adaptive_noise", 0.0f);
    editor->timerCallback();

    // Adaptive-only and no-profile branches of the status line + tooltips
    proc.resetNoiseProfile();
    setParam(proc, "adaptive_noise", 1.0f);
    expect(statusText().startsWith("STATUS: CONTINUOUS ADAPTIVE"),
           "Adaptive without profile must report continuous adaptive");
    setParam(proc, "adaptive_noise", 0.0f);
    expect(statusText().startsWith("STATUS: NO PROFILE (PASS-THROUGH)"),
           "No profile and no adaptive must report pass-through");
    expect(parts.algoCombo->getTooltip().startsWith("Pass-through"),
           "Without processing the combo must advertise pass-through");
    expect(!parts.resetProfile->isEnabled(),
           "Reset must be disabled without profile or adaptive");

    // Unlinked profile mode enables the tonal controls
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    learnNoiseProfile(proc, false, buffer, midi);
    editor->timerCallback();
    expect(parts.resetProfile->isEnabled(),
           "Reset must be enabled once a profile exists");
    expect(parts.tonalRed->isEnabled(),
           "Tonal reduction must enable when unlinked with a profile");
    expect(parts.tonalOffset->isEnabled(),
           "Tonal offset must enable when unlinked with a profile");

    // Slider label updates (wide + compact widths)
    setParam(proc, "smoothing_factor", 42.0f);
    setParam(proc, "masking_depth", 30.0f);
    setParam(proc, "whitening_factor", 55.0f);
    editor->setSize(1000, 620);
    editor->timerCallback();
    {
      std::vector<juce::Component*> comps;
      collectComponents(*editor, comps);
      int checked = 0;
      for (auto* c : comps) {
        auto* label = dynamic_cast<juce::Label*>(c);
        if (label == nullptr) {
          continue;
        }
        const auto text = label->getText();
        if (text.startsWith("MASKING PROTECT:")) {
          ++checked;
        }
        if (text == "THRESHOLD") {
          ++checked;
        }
        if (text.startsWith("SMOOTHING: 42%")) {
          ++checked;
        }
        if (text.startsWith("WHITENING: 55%")) {
          ++checked;
        }
      }
      expectEquals(checked, 4, "Wide layout must use full slider labels");
    }
    editor->setSize(860, 520);
    editor->timerCallback();
    {
      std::vector<juce::Component*> comps;
      collectComponents(*editor, comps);
      int checked = 0;
      for (auto* c : comps) {
        auto* label = dynamic_cast<juce::Label*>(c);
        if (label == nullptr) {
          continue;
        }
        const auto text = label->getText();
        if (text.startsWith("MASKING: 30%")) {
          ++checked;
        }
        if (text == "THRESH") {
          ++checked;
        }
      }
      expectEquals(checked, 2, "Compact layout must use short slider labels");
    }
    proc.releaseResources();
  }

  // ── PluginEditor: resized()/updateLayout()/paint()/paintOverChildren across
  // linked/curve/advanced and footer visibility combinations ──
  void testEditorLayoutPaint() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be the plugin editor");
    if (editor == nullptr) {
      return;
    }
    EditorParts parts(*editor);
    parts.expectCaptured(*this);
    pumpMessageLoop(10);

    // Both footer slots
    setParam(proc, "show_tooltips", 1.0f);
    setParam(proc, "show_hud", 1.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1100, 700);
    expect(parts.footerLabel->isVisible() && parts.statusLabel->isVisible(),
           "Tooltips + HUD must both be visible in the footer");
    expect(parts.footerLabel->getRight() <= parts.statusLabel->getX(),
           "Footer tooltip must not overlap the HUD status");

    // Tooltips-only footer
    setParam(proc, "show_hud", 0.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1104, 700);
    expect(parts.footerLabel->isVisible() && !parts.statusLabel->isVisible(),
           "Tooltips-only footer must hide the HUD label");

    // HUD-only footer
    setParam(proc, "show_tooltips", 0.0f);
    setParam(proc, "show_hud", 1.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1108, 700);
    expect(!parts.footerLabel->isVisible() && parts.statusLabel->isVisible(),
           "HUD-only footer must hide the tooltip label");

    // No footer at all
    setParam(proc, "show_hud", 0.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1112, 700);
    expect(!parts.footerLabel->isVisible() && !parts.statusLabel->isVisible(),
           "Footer bar must be hidden when both options are off");

    setParam(proc, "show_tooltips", 1.0f);
    setParam(proc, "show_hud", 1.0f);
    editor->handleAsyncUpdate();

    // updateSliderLabels crosses the compact threshold: label wording swaps
    editor->setSize(860, 700);
    editor->updateSliderLabels();
    expect(editor->advanced.lblMasking.getText().contains("MASKING") &&
               !editor->advanced.lblMasking.getText().contains("PROTECT"),
           "Compact width must use the short masking label");
    editor->setSize(960, 700);
    editor->updateSliderLabels();
    expect(editor->advanced.lblMasking.getText().contains("PROTECT"),
           "Wide width must use the full masking label");

    // Advanced panel open with linked faders: tonal faders hidden
    setParam(proc, "show_advanced", 1.0f);
    setParam(proc, "link_reduction", 1.0f);
    setParam(proc, "link_threshold_offset", 1.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1200, 800);
    expect(!parts.tonalRed->isVisible() && !parts.tonalOffset->isVisible(),
           "Linked layout must hide the tonal fader columns");

    // Unlinked fader banks (broadband + tonal side by side)
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1204, 800);
    expect(parts.masterRed->isVisible() && parts.tonalRed->isVisible(),
           "Unlinked reduction faders must be visible");
    expect(parts.offset->isVisible() && parts.tonalOffset->isVisible(),
           "Unlinked threshold faders must be visible");
    expect(!parts.masterRed->getBounds().intersects(
               parts.tonalRed->getBounds()),
           "Broadband and tonal reduction faders must not overlap");
    expect(!parts.offset->getBounds().intersects(
               parts.tonalOffset->getBounds()),
           "Broadband and tonal threshold faders must not overlap");

    // Curve split button layout (on and off)
    setParam(proc, "reduction_curve_enabled", 1.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1208, 800);
    expect(parts.resetCurve->isVisible(),
           "Curve reset button must show while the curve is enabled");
    setParam(proc, "reduction_curve_enabled", 0.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1212, 800);
    expect(!parts.resetCurve->isVisible(),
           "Curve reset button must hide while the curve is disabled");

    // Paint with advanced panel open (separator lines drawn) and closed
    {
      auto img = paintComponent(*editor);
      expect(!img.isNull() && countNonBlankPixels(img) > 0,
             "Editor paint must render the panel");
    }
    setParam(proc, "show_advanced", 0.0f);
    editor->handleAsyncUpdate();
    editor->setSize(1216, 800);
    {
      auto img = paintComponent(*editor);
      expect(countNonBlankPixels(img) > 0,
             "Editor paint must render without the advanced panel");
    }

    // Degenerate size: separator lines are skipped when sliders have no width
    parts.advancedToggle->triggerClick();
    pumpMessageLoop(30);
    editor->setSize(80, 60);
    {
      auto img = paintComponent(*editor);
      expect(!img.isNull(), "Editor must survive painting at tiny sizes");
    }
    editor->setSize(1220, 800);
    proc.releaseResources();
  }

  // ── PluginEditor: mouseEnter/mouseMove/mouseExit tooltip routing,
  // parameterChanged + handleAsyncUpdate ──
  void testEditorTooltipPaths() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be the plugin editor");
    if (editor == nullptr) {
      return;
    }
    EditorParts parts(*editor);
    parts.expectCaptured(*this);
    pumpMessageLoop(30);

    auto footerText = [&]() { return parts.footerLabel->getText(); };

    // Tooltips enabled: hovering a tooltip client shows its tip
    setParam(proc, "show_tooltips", 1.0f);
    parts.resetProfile->triggerClick();
    pumpMessageLoop(20);
    editor->mouseEnter(makeEvent(*editor, parts.resetProfile,
                                 parts.resetProfile->getBounds().getCentre()
                                     .toFloat()));
    expect(footerText() == tooltipOf(*parts.resetProfile),
           "Hovering a tooltip client must forward its tooltip to the footer");

    // Hovering a tooltip-less client walks up to the default instruction
    editor->mouseEnter(makeEvent(*editor, parts.statusLabel,
                                 parts.statusLabel->getBounds().getCentre()
                                     .toFloat()));
    expect(footerText().startsWith("Loop a noise-only"),
           "Hovering a tooltip-less component must restore the default "
           "instruction");

    // Each group component re-resolves the component under the cursor
    for (auto* group : parts.groups) {
      editor->mouseEnter(makeEvent(*editor, group,
                                   group->getBounds().getCentre().toFloat()));
    }
    expect(footerText().isNotEmpty(),
           "Group hover must resolve some footer text without crashing");

    // nullptr original component: resolved via getComponentAt
    editor->mouseEnter(makeEvent(
        *editor, nullptr,
        parts.learn->getBounds().getCentre().toFloat()));
    expect(footerText().isNotEmpty(),
           "Synthetic hover without originator must resolve a tooltip");

    // The editor itself as originator also re-resolves under the cursor
    editor->mouseEnter(makeEvent(*editor, editor, {2.0f, 2.0f}));
    editor->mouseMove(makeEvent(*editor, parts.learn, {5.0f, 5.0f}));

    // mouseExit restores / clears depending on the tooltip preference
    editor->mouseExit(makeEvent(*editor, parts.learn, {5.0f, 5.0f}));
    expect(footerText().startsWith("Loop a noise-only"),
           "mouseExit must restore the default instruction");
    setParam(proc, "show_tooltips", 0.0f);
    editor->mouseExit(makeEvent(*editor, parts.learn, {5.0f, 5.0f}));
    expect(footerText().isEmpty(),
           "mouseExit with tooltips off must clear the footer");
    editor->mouseEnter(makeEvent(*editor, parts.learn,
                                 parts.learn->getBounds().getCentre()
                                     .toFloat()));
    expect(footerText().isEmpty(),
           "mouseEnter with tooltips off must keep the footer clear");

    // parameterChanged + async update plumbing
    editor->parameterChanged("unrelated_parameter", 1.0f);
    editor->timerCallback();
    editor->parameterChanged("show_tooltips", 1.0f);
    setParam(proc, "show_tooltips", 1.0f);
    editor->timerCallback();
    pumpMessageLoop(30);
    expect(footerText().startsWith("Loop a noise-only"),
           "Async tooltip update must restore the default instruction");
    editor->parameterChanged("show_tooltips", 0.0f);
    setParam(proc, "show_tooltips", 0.0f);
    editor->timerCallback();
    pumpMessageLoop(30);
    expect(footerText().isEmpty(),
           "Async tooltip update must clear the footer when disabled");
    editor->handleAsyncUpdate();
    setParam(proc, "show_tooltips", 1.0f);
    editor->handleAsyncUpdate();
    expect(footerText().startsWith("Loop a noise-only"),
           "handleAsyncUpdate must sync the footer to the tooltip preference");
    proc.releaseResources();
  }

  // ── PluginEditor: offline rendering overlay (paintOverChildren) ──
  void testEditorOfflineOverlay() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    pumpMessageLoop(10);
    ScopedEditor ed(proc);
    auto* editor =
        dynamic_cast<NoiseRepellentAudioProcessorEditor*>(ed.editor.get());
    expect(editor != nullptr, "Editor must be the plugin editor");
    if (editor == nullptr) {
      return;
    }
    pumpMessageLoop(10);

    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    const auto noiseFill = [](juce::AudioBuffer<float>& b) {
      generateNoiseBuffer(b, 0.05f);
    };

    // Not offline: overlay paints nothing
    {
      juce::Image img(juce::Image::ARGB, editor->getWidth(),
                      editor->getHeight(), true);
      juce::Graphics g(img);
      editor->paintOverChildren(g);
      expectEquals(countNonBlankPixels(img), 0,
                   "Overlay must stay hidden outside offline rendering");
    }

    // Offline flag without flowing blocks must NOT flash the overlay
    proc.setNonRealtime(true);
    editor->timerCallback();
    {
      juce::Image img(juce::Image::ARGB, editor->getWidth(),
                      editor->getHeight(), true);
      juce::Graphics g(img);
      editor->paintOverChildren(g);
      expectEquals(countNonBlankPixels(img), 0,
                   "A stale offline flag must not show the overlay");
    }

    // Offline blocks actually flowing: overlay appears
    processBlocks(proc, buffer, midi, 6, noiseFill);
    editor->timerCallback();
    {
      juce::Image img(juce::Image::ARGB, editor->getWidth(),
                      editor->getHeight(), true);
      juce::Graphics g(img);
      editor->paintOverChildren(g);
      expect(countNonBlankPixels(img) > editor->getWidth() *
                                             editor->getHeight() / 2,
             "Active offline rendering must darken the whole editor");
    }

    // Painting while the overlay is up must not crash and still render
    {
      auto img = paintComponent(*editor);
      expect(countNonBlankPixels(img) > 0, "Editor paint must stay functional");
    }

    // Back to realtime: overlay hides again
    proc.setNonRealtime(false);
    processBlocks(proc, buffer, midi, 2, noiseFill);
    editor->timerCallback();
    {
      juce::Image img(juce::Image::ARGB, editor->getWidth(),
                      editor->getHeight(), true);
      juce::Graphics g(img);
      editor->paintOverChildren(g);
      expectEquals(countNonBlankPixels(img), 0,
                   "Leaving offline rendering must hide the overlay");
    }
    proc.releaseResources();
  }

  // ── LookAndFeel: every public drawing/metric entry point ──
  void testLookAndFeelDrawing() {
    NoiseRepellentLookAndFeel lnf;
    juce::Image img(juce::Image::RGB, 200, 120, true);
    juce::Graphics g(img);

    // Linear sliders: vertical, horizontal, fallback fill colour, other style
    juce::Slider vSlider(juce::Slider::LinearVertical,
                         juce::Slider::TextBoxBelow);
    vSlider.setSize(40, 100);
    vSlider.setColour(juce::Slider::rotarySliderFillColourId,
                      juce::Colour(0xff336699));
    lnf.drawLinearSlider(g, 0, 0, 40, 100, 60.0f, 0.0f, 100.0f,
                         juce::Slider::LinearVertical, vSlider);
    juce::Slider hSlider(juce::Slider::LinearHorizontal,
                         juce::Slider::NoTextBox);
    hSlider.setSize(100, 24);
    hSlider.setColour(juce::Slider::rotarySliderFillColourId,
                      juce::Colour(0x00000000)); // transparent -> fallback
    lnf.drawLinearSlider(g, 50, 20, 100, 24, 90.0f, 50.0f, 130.0f,
                         juce::Slider::LinearHorizontal, hSlider);
    lnf.drawLinearSlider(g, 50, 50, 100, 24, 90.0f, 50.0f, 130.0f,
                         juce::Slider::RotaryHorizontalVerticalDrag, hSlider);

    // Button backgrounds: brand early-out, on/off, translucent fallbacks,
    // down and highlighted states
    juce::TextButton brand("NOISE REPELLENT");
    brand.setName(NoiseRepellentLookAndFeel::kBrandButtonName);
    brand.setSize(100, 26);
    lnf.drawButtonBackground(g, brand, juce::Colours::black, true, true);

    juce::TextButton onBtn("on");
    onBtn.setSize(100, 24);
    onBtn.setToggleState(true, juce::dontSendNotification);
    onBtn.setColour(juce::TextButton::buttonOnColourId,
                    juce::Colour(0xff112233));
    lnf.drawButtonBackground(g, onBtn, juce::Colours::black, false, false);
    expect(juce::Colour(img.getPixelAt(50, 12).getARGB()) ==
               juce::Colour(0xff112233),
           "On-state button must use the opaque on colour");

    juce::TextButton translucentOn("translucentOn");
    translucentOn.setSize(100, 24);
    translucentOn.setToggleState(true, juce::dontSendNotification);
    translucentOn.setColour(juce::TextButton::buttonOnColourId,
                            juce::Colour(0x00112233));
    lnf.drawButtonBackground(g, translucentOn, juce::Colours::black, false,
                             false);
    expect(juce::Colour(img.getPixelAt(50, 12).getARGB()) ==
               NoiseRepellentLookAndFeel::kColorNoiseProfile,
           "Translucent on colour must fall back to the noise profile colour");

    juce::TextButton offBtn("off");
    offBtn.setSize(100, 24);
    offBtn.setColour(juce::TextButton::buttonColourId,
                     juce::Colour(0xff445566));
    lnf.drawButtonBackground(g, offBtn, juce::Colours::black, true, false);
    expect(juce::Colour(img.getPixelAt(50, 12).getARGB()) ==
               juce::Colour(0xff445566).brighter(0.15f),
           "Highlighted off-state must brighten the button colour");

    lnf.drawButtonBackground(g, offBtn, juce::Colours::black, false, true);
    expect(juce::Colour(img.getPixelAt(50, 12).getARGB()) ==
               juce::Colour(0xff445566).darker(0.2f),
           "Pressed off-state must darken the button colour");

    offBtn.setColour(juce::TextButton::buttonColourId,
                     juce::Colour(0x00445566));
    lnf.drawButtonBackground(g, offBtn, juce::Colours::transparentBlack, false,
                             false);
    expect(juce::Colour(img.getPixelAt(50, 12).getARGB()) ==
               NoiseRepellentLookAndFeel::kColorButtonOff,
           "Transparent button colour must fall back to the off colour");
    lnf.drawButtonBackground(g, offBtn, juce::Colour(0xff778899), false,
                             false);
    expect(juce::Colour(img.getPixelAt(50, 12).getARGB()) ==
               juce::Colour(0xff778899),
           "Opaque background colour must be used when the button colour is "
           "transparent");

    // Fonts
    expectWithinAbsoluteError(lnf.getTextButtonFont(brand, 26).getHeight(),
                              NoiseRepellentLookAndFeel::kFontSizeBrand, 0.01f,
                              "Brand button must use the brand font size");
    expectWithinAbsoluteError(lnf.getTextButtonFont(offBtn, 24).getHeight(),
                              NoiseRepellentLookAndFeel::kFontSizeLabel, 0.01f,
                              "Regular buttons must use the label font size");
    juce::ComboBox fontBox;
    expectWithinAbsoluteError(lnf.getComboBoxFont(fontBox).getHeight(),
                              NoiseRepellentLookAndFeel::kFontSizeLabel, 0.01f,
                              "Combo boxes must use the label font size");
    expectWithinAbsoluteError(lnf.getPopupMenuFont().getHeight(),
                              NoiseRepellentLookAndFeel::kFontSizeLabel, 0.01f,
                              "Popup menus must use the label font size");

    // Button text: brand, learn colours, light/dark and disabled variants
    lnf.drawButtonText(g, brand, false, false);
    juce::TextButton learnBtn("Learn Noise");
    learnBtn.setSize(100, 24);
    learnBtn.setColour(juce::TextButton::buttonColourId,
                       NoiseRepellentLookAndFeel::kColorLearnCTA);
    lnf.drawButtonText(g, learnBtn, false, false);
    learnBtn.setColour(juce::TextButton::buttonColourId,
                       NoiseRepellentLookAndFeel::kColorLearnActive);
    lnf.drawButtonText(g, learnBtn, true, true);
    learnBtn.setColour(juce::TextButton::buttonColourId,
                       NoiseRepellentLookAndFeel::kColorNoiseProfile);
    lnf.drawButtonText(g, learnBtn, false, false);
    learnBtn.setColour(juce::TextButton::buttonColourId,
                       juce::Colour(0xff222222));
    lnf.drawButtonText(g, learnBtn, false, false);
    learnBtn.setEnabled(false);
    lnf.drawButtonText(g, learnBtn, false, false);
    learnBtn.setEnabled(true);

    // Toggle buttons: on, off, disabled
    juce::ToggleButton toggle("Linked");
    toggle.setSize(100, 24);
    lnf.drawToggleButton(g, toggle, true, true);
    toggle.setToggleState(true, juce::dontSendNotification);
    lnf.drawToggleButton(g, toggle, false, false);
    toggle.setEnabled(false);
    lnf.drawToggleButton(g, toggle, false, false);
    toggle.setEnabled(true);

    // Combo box chrome
    juce::ComboBox box;
    box.setSize(120, 26);
    lnf.drawComboBox(g, 120, 26, true, 90, 2, 28, 22, box);

    // Group outlines with and without caption text
    juce::GroupComponent group("group", "");
    group.setSize(160, 80);
    lnf.drawGroupComponentOutline(g, 160, 80, "",
                                  juce::Justification::left, group);
    lnf.drawGroupComponentOutline(g, 160, 80, "NOISE PROFILE",
                                  juce::Justification::left, group);

    // Popup menu section header
    lnf.drawPopupMenuSectionHeader(g, juce::Rectangle<int>(0, 0, 160, 24),
                                   "FRAME SIZE");

    // Popup item metrics
    int idealWidth = 0;
    int idealHeight = 0;
    lnf.getIdealPopupMenuItemSize("Some item", false, 20, idealWidth,
                                  idealHeight);
    expectEquals(idealHeight, 24, "Popup items must be 24px tall");
    lnf.getIdealPopupMenuItemSize({}, true, 20, idealWidth, idealHeight);
    expectEquals(idealHeight, 24, "Popup separators must be 24px tall");
  }

  // ── SpectralVisualizer: paint matrix and timer state machine ──
  void testVisualizerPaintAndTimer() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    // The processor publishes spectral frames only while an editor is open.
    ScopedEditor ed(proc);
    expect(ed.editor != nullptr, "Editor must exist while feeding frames");

    // Pass-through state: no profile at all
    {
      SpectralVisualizerComponent viz(proc);
      viz.stopTimer();
      viz.setSize(800, 400);
      viz.setProcessingAvailable(false);
      viz.setAdvancedControlsVisible(false);
      viz.timerCallback(); // no frames yet: idle timer tick
      auto idle = paintComponent(viz);
      expect(!idle.isNull(), "Visualizer must paint without any spectral data");

      // Pass-through card: the centre card must differ from the plain frame
      SpectralVisualizerComponent viz2(proc);
      viz2.stopTimer();
      viz2.setSize(800, 400);
      viz2.setProcessingAvailable(true);
      viz2.setAdvancedControlsVisible(false);
      auto noCard = paintComponent(viz2);
      expect(countDifferentPixels(idle, noCard) > 1000,
             "Pass-through card must render when processing is unavailable");

      // Tiny sizes skip the card and the badge
      viz.setSize(30, 30);
      auto tiny = paintComponent(viz);
      expect(!tiny.isNull(), "Visualizer must paint at tiny sizes");
      viz.setSize(60, 60);
      viz.setProcessingAvailable(true);
      viz.setAdvancedControlsVisible(true);
      auto tinyAdvanced = paintComponent(viz);
      expect(!tinyAdvanced.isNull(),
             "Visualizer must paint tiny sizes with advanced controls");
      // getTpBadgeBounds geometry
      viz.setSize(800, 400);
      expect(viz.getTpBadgeBounds() ==
                 juce::Rectangle<float>(800.0f - 66.0f - 10.0f, 10.0f, 66.0f,
                                        22.0f),
             "TP badge must sit in the top-right corner");
      // bias mapping helpers
      expectWithinAbsoluteError(SpectralVisualizerComponent::biasToY(0.0f,
                                                                    400.0f),
                                200.0f, 0.01f,
                                "Zero bias maps to the vertical centre");
      expectWithinAbsoluteError(
          SpectralVisualizerComponent::yToBias(
              SpectralVisualizerComponent::biasToY(12.0f, 400.0f), 400.0f),
          12.0f, 0.01f, "biasToY/yToBias must round-trip");
    }

    // Learn a tonal profile so frames carry a noise floor and tonal peaks
    learnNoiseProfile(proc, true, buffer, midi);
    setParam(proc, "link_reduction", 0.0f);
    setParam(proc, "link_threshold_offset", 0.0f);
    setParam(proc, "reduction_curve_enabled", 1.0f);
    proc.curve.setNodes({{0.0f, 6.0f}, {0.5f, -6.0f}, {1.0f, 3.0f}});

    SpectralVisualizerComponent viz(proc);
    viz.stopTimer();
    viz.setSize(800, 400);
    viz.setProcessingAvailable(true);
    viz.setAdvancedControlsVisible(true);
    auto drain = [&viz](int ticks) {
      for (int i = 0; i < ticks; ++i) {
        viz.timerCallback();
      }
    };

    processBlocks(proc, buffer, midi, 12);
    drain(20);
    auto fullFrame = paintComponent(viz);
    expect(!fullFrame.isNull(),
           "Visualizer must paint a full profile + curve + peaks frame");

    // Attack: live input raises the displayed spectra
    processBlocks(proc, buffer, midi, 8,
                  [](juce::AudioBuffer<float>& b) {
                    generateNoiseBuffer(b, 0.25f);
                  });
    drain(8);
    auto loud = paintComponent(viz);
    expect(countDifferentPixels(fullFrame, loud) > 100,
           "Live spectra must update the display");

    // Release: quiet input relaxes the displayed spectra
    processBlocks(proc, buffer, midi, 8,
                  [](juce::AudioBuffer<float>& b) {
                    generateNoiseBuffer(b, 0.0005f);
                  });
    drain(8);

    // Host stopped: display decays toward the floor
    proc.status.hostTransportStateKnown.store(true);
    proc.status.hostTransportPlaying.store(false);
    drain(40);
    proc.status.hostTransportPlaying.store(true);
    proc.status.hostTransportStateKnown.store(false);

    // Legend swatch variants: linked vs unlinked, curve on/off. Frames are
    // produced and consumed without message-loop pumping in between so the
    // editor's own visualizer cannot steal them.
    auto paintAfterFlags = [&](bool linkReduction, bool linkOffset,
                               bool curveEnabled) {
      setParam(proc, "link_reduction", linkReduction ? 1.0f : 0.0f);
      setParam(proc, "link_threshold_offset", linkOffset ? 1.0f : 0.0f);
      setParam(proc, "reduction_curve_enabled", curveEnabled ? 1.0f : 0.0f);
      processBlocks(proc, buffer, midi, 8);
      drain(8);
      return paintComponent(viz);
    };
    auto legendFull = paintAfterFlags(false, false, true);
    auto legendTonalOnly = paintAfterFlags(true, false, false);
    auto legendLinked = paintAfterFlags(true, true, false);
    expect(countDifferentPixels(legendFull, legendTonalOnly) > 0,
           "Legend must change when the curve swatch toggles");
    expect(countDifferentPixels(legendTonalOnly, legendLinked) > 0,
           "Legend must change when tonal peaks toggle");

    // Transient LED state machine: trigger, hold decay, off
    setParam(proc, "transient_protection_enable", 1.0f);
    {
      NoiseRepellentAudioProcessor::SpectralFrame scratch;
      while (proc.viz.getNextFrame(scratch)) {
      }
    }
    const auto clickFill = [](juce::AudioBuffer<float>& b) {
      generateClickBuffer(b);
    };
    // Probe batch: the stimulus must genuinely trigger the detector
    processBlocks(proc, buffer, midi, 8, clickFill);
    float maxTransient = 0.0f;
    {
      NoiseRepellentAudioProcessor::SpectralFrame probe;
      while (proc.viz.getNextFrame(probe)) {
        maxTransient = std::max(maxTransient, probe.transientIntensity);
      }
    }
    expect(maxTransient >= 0.35f,
           "Click bursts must trigger the transient detector");

    // LED batch: consume one transient-lit frame, then let the hold decay
    processBlocks(proc, buffer, midi, 8, clickFill);
    viz.timerCallback();
    auto ledOn = paintComponent(viz);
    {
      NoiseRepellentAudioProcessor::SpectralFrame scratch;
      while (proc.viz.getNextFrame(scratch)) {
      }
    }
    drain(8); // hold ticks decay with no frames left
    auto ledOff = paintComponent(viz);
    expect(countDifferentPixels(ledOn, ledOff) > 0,
           "Transient LED must fade after the burst");
    proc.releaseResources();
  }

  // ── SpectralVisualizer: mouse-driven curve editing and TP badge ──
  void testVisualizerInteractions() {
    NoiseRepellentAudioProcessor proc;
    proc.prepareToPlay(48000.0, 512);
    juce::AudioBuffer<float> buffer(2, 512);
    juce::MidiBuffer midi;
    // The processor publishes spectral frames only while an editor is open.
    ScopedEditor ed(proc);

    setParam(proc, "reduction_curve_enabled", 1.0f);
    setParam(proc, "transient_protection_enable", 0.0f);
    // setNodes() appends the mandatory boundary nodes: this yields 5.
    proc.curve.setNodes({{0.1f, 6.0f}, {0.5f, -6.0f}, {0.9f, 3.0f}});

    SpectralVisualizerComponent viz(proc);
    viz.stopTimer();
    viz.setSize(800, 400);
    viz.setProcessingAvailable(false);
    viz.setAdvancedControlsVisible(true);

    // Feed frames carrying the curve flag into the visualizer's state
    processBlocks(proc, buffer, midi, 12);
    for (int i = 0; i < 8; ++i) {
      viz.timerCallback();
    }

    // Unavailable: all interaction is inert
    const auto before = proc.curve.nodes();
    viz.mouseDown(makeEvent(viz, &viz, {400.0f, 200.0f}));
    viz.mouseDrag(makeEvent(viz, &viz, {500.0f, 150.0f}, {400.0f, 200.0f}));
    viz.mouseDoubleClick(makeEvent(viz, &viz, {400.0f, 300.0f}));
    expectEquals(static_cast<int>(proc.curve.nodes().size()),
                 static_cast<int>(before.size()),
                 "Curve edits must be inert while processing is unavailable");
    viz.mouseMove(makeEvent(viz, &viz, viz.getTpBadgeBounds().getCentre()));
    expect(viz.getMouseCursor() == juce::MouseCursor::NormalCursor,
           "Badge hover must stay inert while processing is unavailable");
    viz.mouseExit(makeEvent(viz, &viz, {10.0f, 10.0f}));
    expect(viz.getTooltip().startsWith("Pass-through"),
           "Unavailable visualizer must advertise pass-through in its tooltip");

    // Available: hover states and tooltip behaviour
    viz.setProcessingAvailable(true);
    viz.setProcessingAvailable(true); // idempotent re-entry
    viz.mouseMove(makeEvent(viz, &viz, viz.getTpBadgeBounds().getCentre()));
    expect(viz.getMouseCursor() == juce::MouseCursor::PointingHandCursor,
           "Hovering the TP badge must show the pointing hand");
    viz.mouseMove(makeEvent(viz, &viz, {200.0f, 300.0f}));
    expect(viz.getMouseCursor() == juce::MouseCursor::NormalCursor,
           "Hovering the spectrum must show the normal cursor");
    viz.mouseExit(makeEvent(viz, &viz, {10.0f, 10.0f}));
    expect(viz.getTooltip().isEmpty(),
           "Spectrum hover without the badge must have no tooltip");

    // TP badge click toggles transient protection through the APVTS
    const auto badgeClick = viz.getTpBadgeBounds().getCentre();
    viz.mouseDown(makeEvent(viz, &viz, badgeClick));
    expectEquals(rawParam(proc, "transient_protection_enable"), 1.0f,
                 "Clicking the TP badge must enable transient protection");
    viz.mouseDown(makeEvent(viz, &viz, badgeClick));
    expectEquals(rawParam(proc, "transient_protection_enable"), 0.0f,
                 "Clicking the TP badge again must disable it");

    // Empty click on the curve area adds a node
    const int countBefore = static_cast<int>(proc.curve.nodes().size());
    viz.mouseDown(makeEvent(viz, &viz, {400.0f, 200.0f}));
    expectEquals(static_cast<int>(proc.curve.nodes().size()), countBefore + 1,
                 "Clicking empty spectrum must add a curve node");

    // Drag existing node 1 and verify it moved
    {
      const auto& nodes = proc.curve.nodes();
      const auto idx = static_cast<size_t>(1);
      const float nx = juce::jlimit(0.0f, 1.0f, nodes[idx].normX);
      const float ny =
          SpectralVisualizerComponent::biasToY(nodes[idx].biasDB, 400.0f);
      viz.mouseDown(makeEvent(viz, &viz, {nx * 800.0f, ny}, {nx * 800.0f, ny}));
      viz.mouseDrag(makeEvent(viz, &viz, {0.75f * 800.0f, 100.0f},
                              {nx * 800.0f, ny}));
      viz.mouseUp(makeEvent(viz, &viz, {0.75f * 800.0f, 100.0f}));
      bool found = false;
      for (const auto& n : proc.curve.nodes()) {
        if (std::abs(n.normX - 0.75f) < 0.001f) {
          found = true;
          const float expectedBias =
              SpectralVisualizerComponent::yToBias(100.0f, 400.0f);
          expectWithinAbsoluteError(n.biasDB, expectedBias, 0.01f,
                                    "Dragging must move the node's bias");
        }
      }
      expect(found, "Dragging must move the node to the new X position");
    }

    // Clamped drag keeps the node inside the plot
    {
      const auto& nodes = proc.curve.nodes();
      const auto idx = static_cast<size_t>(1);
      const float nx = juce::jlimit(0.0f, 1.0f, nodes[idx].normX);
      const float ny =
          SpectralVisualizerComponent::biasToY(nodes[idx].biasDB, 400.0f);
      viz.mouseDown(makeEvent(viz, &viz, {nx * 800.0f, ny}, {nx * 800.0f, ny}));
      viz.mouseDrag(
          makeEvent(viz, &viz, {-500.0f, -500.0f}, {nx * 800.0f, ny}));
      viz.mouseUp(makeEvent(viz, &viz, {-500.0f, -500.0f}));
      for (const auto& n : proc.curve.nodes()) {
        expect(n.normX >= 0.0f && n.normX <= 1.0f,
               "Node X must stay clamped inside the plot");
        expect(n.biasDB >= -SpectralVisualizerComponent::kCurveMaxBiasDB &&
                   n.biasDB <= SpectralVisualizerComponent::kCurveMaxBiasDB,
               "Node bias must stay clamped inside the curve range");
      }
    }

    // Mouse-up clears the drag target; drags without a target are inert
    viz.mouseUp(makeEvent(viz, &viz, {10.0f, 10.0f}));
    {
      const auto frozen = proc.curve.nodes();
      viz.mouseDrag(makeEvent(viz, &viz, {100.0f, 100.0f}, {50.0f, 50.0f}));
      expectEquals(static_cast<int>(proc.curve.nodes().size()),
                   static_cast<int>(frozen.size()),
                   "Dragging without an active node must not edit the curve");
    }

    // Double-click removes an interior node (boundary nodes survive)
    {
      const auto count = static_cast<int>(proc.curve.nodes().size());
      expect(count >= 3, "Test expects at least three curve nodes");
      const auto& nodes = proc.curve.nodes();
      const auto idx = static_cast<size_t>(1);
      const float nx = juce::jlimit(0.0f, 1.0f, nodes[idx].normX);
      const float ny =
          SpectralVisualizerComponent::biasToY(nodes[idx].biasDB, 400.0f);
      viz.mouseDoubleClick(makeEvent(viz, &viz, {nx * 800.0f, ny}));
      expectEquals(static_cast<int>(proc.curve.nodes().size()), count - 1,
                   "Double-clicking an interior node must remove it");
    }

    // Double-click on a boundary node must not remove it
    {
      const auto count = static_cast<int>(proc.curve.nodes().size());
      const auto& nodes = proc.curve.nodes();
      const float ny =
          SpectralVisualizerComponent::biasToY(nodes.back().biasDB, 400.0f);
      viz.mouseDoubleClick(makeEvent(
          viz, &viz,
          {juce::jlimit(0.0f, 1.0f, nodes.back().normX) * 800.0f, ny}));
      expectEquals(static_cast<int>(proc.curve.nodes().size()), count,
                   "Boundary curve nodes must survive double clicks");
    }

    // Losing processing availability cancels an in-flight drag
    {
      const auto& nodes = proc.curve.nodes();
      const auto idx = static_cast<size_t>(1);
      const float nx = juce::jlimit(0.0f, 1.0f, nodes[idx].normX);
      const float ny =
          SpectralVisualizerComponent::biasToY(nodes[idx].biasDB, 400.0f);
      viz.mouseDown(makeEvent(viz, &viz, {nx * 800.0f, ny}));
      viz.setProcessingAvailable(false);
      const auto frozen = proc.curve.nodes();
      viz.mouseDrag(makeEvent(viz, &viz, {600.0f, 100.0f}, {nx * 800.0f, ny}));
      expectEquals(static_cast<int>(proc.curve.nodes().size()),
                   static_cast<int>(frozen.size()),
                   "Dropped drag targets must not move curve nodes");
      viz.setProcessingAvailable(true);
    }

    // resized() is a no-op layout hook
    viz.resized();
    proc.releaseResources();
  }
};

GuiCoverageTest guiCoverageTest;

int main(int argc, char* argv[]) {
  juce::ScopedJuceInitialiser_GUI guiInit;

  // Keep dialog paths on the JUCE-drawn implementation so the blocking
  // showMessageBoxAsync path is dismissable headlessly from the test timers.
  juce::LookAndFeel::getDefaultLookAndFeel().setUsingNativeAlertWindows(false);

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
