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

#pragma once

#include "GUI/LookAndFeel.h"
#include "GUI/SpectralVisualizer.h"
#include "PluginProcessor.h"
#include <atomic>
#include <functional>
#include <juce_gui_basics/juce_gui_basics.h>

// Headless-testable content builders: they assemble menu items / dialog text
// without showing native windows (which crash under Linux CI Xvfb).
juce::PopupMenu buildEditorPreferencesMenu(
    NoiseRepellentAudioProcessor& processor,
    juce::LookAndFeel* menuLookAndFeel,
    const std::function<void()>& relayout);
juce::PopupMenu buildAdaptiveMethodMenu(int currentMethod,
                                        juce::LookAndFeel* menuLookAndFeel);
juce::String buildAboutMessage();

class NoiseRepellentAudioProcessorEditor
    : public juce::AudioProcessorEditor,
      public juce::Timer,
      public juce::AudioProcessorValueTreeState::Listener,
      public juce::AsyncUpdater {
public:
  explicit NoiseRepellentAudioProcessorEditor(NoiseRepellentAudioProcessor&);
  ~NoiseRepellentAudioProcessorEditor() override;

  void paint(juce::Graphics&) override;
  void paintOverChildren(juce::Graphics&) override;
  void resized() final;
  void timerCallback() override;
  void mouseEnter(const juce::MouseEvent&) override;
  void mouseMove(const juce::MouseEvent&) override;
  void mouseExit(const juce::MouseEvent&) override;

  // Applies an estimation-method choice made in the popup menu
  // (headless-testable; showing the menu itself needs a display).
  void applyAdaptiveMethodChoice(int result);

  void parameterChanged(const juce::String& parameterID,
                        float newValue) override;
  void handleAsyncUpdate() override;

private:
  // Parameter Attachments
  using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
  using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;
  using ComboBoxAttachment =
      juce::AudioProcessorValueTreeState::ComboBoxAttachment;

  friend class GuiCoverageTest; // headless coverage tests poke internal controls

  // Controls are grouped in one struct per UI section
  struct HeaderControls {
    juce::TextButton brandButton;
    juce::TextButton btnPreferences{juce::CharPointer_UTF8("\xe2\x96\xbc")};
    juce::Label lblAlgoHeader{"lblAlgoHeader", "SMOOTHING METHOD"};
    juce::ComboBox comboAlgoMode;
    juce::TextButton btnAdvancedToggle{"ADVANCED"};
    juce::ToggleButton btnDelta{"Delta"};
    juce::ToggleButton btnBypass{"Bypass"};
  };

  struct ProfileControls {
    juce::GroupComponent groupProfile{"groupProfile", "NOISE PROFILE"};
    juce::TextButton btnLearn{"Learn Noise"};
    juce::TextButton btnAdaptiveNoise{"Adaptive"};
    juce::TextButton btnAdaptiveArrow{juce::CharPointer_UTF8("\xe2\x96\xbc")};
    juce::TextButton btnResetProfile{juce::CharPointer_UTF8("\xe2\x86\xba")};
    juce::Slider sliderOffset{juce::Slider::LinearVertical,
                              juce::Slider::TextBoxBelow};
    juce::Label lblOffset{"lblOffset", "THRESHOLD"};
    juce::ToggleButton btnLinkOffset{"Linked"};
    juce::Label lblMasterOffset{"lblMasterOffset", "BROADBAND"};
    juce::Label lblTonalOffset{"lblTonalOffset", "TONAL"};
    juce::Slider sliderTonalOffset{juce::Slider::LinearVertical,
                                   juce::Slider::TextBoxBelow};
    juce::Label lblProfileStatus;
  };

  struct DenoisingControls {
    juce::GroupComponent groupDenoising{"groupDenoising",
                                        "DENOISING PROCESSING"};
    juce::Label lblReductionHeader{"lblReductionHeader", "REDUCTION"};
    juce::ToggleButton btnLink{"Linked"};
    juce::ToggleButton btnCurveToggle{"Curve"};
    juce::TextButton btnResetCurve{juce::CharPointer_UTF8("\xe2\x86\xba")};
    juce::Label lblMasterRed{"lblMasterRed", "REDUCTION"};
    juce::Label lblTonalRed{"lblTonalRed", "TONAL"};
    juce::Slider sliderMasterRed{juce::Slider::LinearVertical,
                                 juce::Slider::TextBoxBelow};
    juce::Slider sliderTonalRed{juce::Slider::LinearVertical,
                                juce::Slider::TextBoxBelow};
  };

  struct AdvancedControls {
    juce::GroupComponent groupAdvanced{"groupAdvanced", "ADVANCED CONTROLS"};
    juce::ComboBox comboMethod;
    juce::Slider sliderSmoothing{juce::Slider::LinearHorizontal,
                                 juce::Slider::NoTextBox};
    juce::Slider sliderMasking{juce::Slider::LinearHorizontal,
                               juce::Slider::NoTextBox};
    juce::Slider sliderWhitening{juce::Slider::LinearHorizontal,
                                 juce::Slider::NoTextBox};
    juce::Slider sliderAggressiveness{juce::Slider::LinearHorizontal,
                                      juce::Slider::NoTextBox};
    juce::Label lblMethod{"lblMethod", "ESTIMATION METHOD"};
    juce::Label lblSmoothing{"lblSmoothing", "SMOOTHING"};
    juce::Label lblMasking{"lblMasking", "MASKING PROTECT"};
    juce::Label lblWhitening{"lblWhitening", "WHITENING"};
    juce::Label lblAggressiveness{"lblAggressiveness", "AGGRESSIVENESS"};
  };

  struct Attachments {
    std::unique_ptr<ComboBoxAttachment> attachAlgoMode;
    std::unique_ptr<ButtonAttachment> attachLearn;
    std::unique_ptr<ButtonAttachment> attachAdaptive;
    std::unique_ptr<ButtonAttachment> attachLink;
    std::unique_ptr<ButtonAttachment> attachLinkOffset;
    std::unique_ptr<ButtonAttachment> attachCurveToggle;
    std::unique_ptr<ButtonAttachment> attachDelta;
    std::unique_ptr<ButtonAttachment> attachBypass;
    std::unique_ptr<ButtonAttachment> attachShowAdvanced;
    std::unique_ptr<ComboBoxAttachment> attachMethod;
    std::unique_ptr<SliderAttachment> attachMasterRed;
    std::unique_ptr<SliderAttachment> attachTonalRed;
    std::unique_ptr<SliderAttachment> attachOffset;
    std::unique_ptr<SliderAttachment> attachTonalOffset;
    std::unique_ptr<SliderAttachment> attachSmoothing;
    std::unique_ptr<SliderAttachment> attachMasking;
    std::unique_ptr<SliderAttachment> attachWhitening;
    std::unique_ptr<SliderAttachment> attachAggressiveness;
  };

  NoiseRepellentAudioProcessor& audioProcessor;
  NoiseRepellentLookAndFeel customLookAndFeel;

  std::atomic<bool> tooltipStateDirty{false};

  HeaderControls header;
  ProfileControls profile;
  DenoisingControls denoising;
  AdvancedControls advanced;
  Attachments attachments;

  SpectralVisualizerComponent spectralVisualizer;

  // Footer Tooltip Bar
  juce::Label footerTooltipLabel;

  bool isAdvancedVisible = true;
  bool wasOfflineRendering = false;
  // Last observed offline block count — the overlay shows only while this
  // advances (real render in progress), not on a stale host offline flag.
  uint64_t lastSeenOfflineBlockCount = 0;

  void updateLayout();
  void updateSliderLabels();
  void updateProfileStatus();
  void showAboutBox() const;

  // Constructor-time setup split out of the constructor body
  void showPreferencesMenu();
  void showAdaptiveMethodMenu();
  void handleResetProfileClick();
  void handleLinkToggled();
  void handleLinkOffsetToggled();

  // updateProfileStatus stages
  void applyHeaderStates(bool controlsActive, bool pluginActive,
                         bool processingAvailable, bool lowLatency);
  void applyProfileStates(bool pluginActive, bool controlsActive,
                          bool processingAvailable, bool isAdaptive,
                          bool hasProfile);
  void applyReductionStates(bool controlsActive, bool processingAvailable,
                            bool lowLatency, bool isAdaptive, bool hasProfile);
  void applyOffsetStates(bool controlsActive, bool processingAvailable,
                         bool lowLatency, bool isAdaptive, bool hasProfile);
  void applyAdvancedStates(bool controlsActive, bool processingAvailable);

  void layoutFooterBar(juce::Rectangle<int>& area);

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(
      NoiseRepellentAudioProcessorEditor)
};
