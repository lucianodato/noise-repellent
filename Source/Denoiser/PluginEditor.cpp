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

#include "PluginEditor.h"
#include "PluginProcessor.h"
#include "specbleach_version.h"
#include <array>

namespace {

// Shared tooltip/instruction strings (single source; previously copy-pasted)
const juce::String kTipDefaultInstruction =
    "Loop a noise-only segment in your DAW and click Learn Noise\nto "
    "capture a profile, then adjust reduction level.";
const juce::String kTipAwaitingProfile =
    "Pass-through: learn a noise profile or enable Adaptive to activate "
    "denoising controls.";
const juce::String kTipAdaptive =
    "Toggle continuous adaptive noise estimation.\nWorks standalone or on "
    "top of a manually learned profile.";
const juce::String kTipAggressiveness =
    "Adjust noise profile aggressiveness using statistical variance.\n"
    "(-1.0: Median, 0.0: Mean, +1.0: Mean + 2 Standard Deviations).";
const juce::String kTipLink =
    "Link broadband and tonal reduction controls together\nfor unified "
    "adjustment.";
const juce::String kTipResetCurve = "Reset reduction curve to flat 0 dB line.";
const juce::String kTipThreshold =
    "Shift noise profile threshold up or down\nin decibels (-12 to +12 dB).";
// Header
const juce::String kTipAbout = "About Noise Repellent";
const juce::String kTipPreferences = "Preferences";
const juce::String kTipPreferencesMenu = "Plugin preferences menu.";
const juce::String kTipAlgoMode = juce::String::fromUTF8(
    "Four methods, same controls — none outranks the others.\n"
    "Temporal: fast; NLM: detail; + DFTT: tonal; BM3D: steady.");
const juce::String kTipAdvancedToggle =
    "Toggle Advanced DSP Controls (Smoothing, Masking, Whitening, "
    "Aggressiveness).";
const juce::String kTipAdvancedShow =
    "Show Advanced DSP Controls (Smoothing, Masking Protect, Whitening,\n "
    "Bias Curve, Tonal Split & Aggressiveness).";
// Noise profile box
const juce::String kTipLearnInitial =
    "Capture static noise profile from current audio input\n(supports "
    "multi-section accumulation).";
const juce::String kTipLearnCapturing =
    "Capturing noise profile from current audio playback...";
const juce::String kTipLearnAccumulate =
    "Capture and accumulate an additional noise section\ninto the current "
    "profile.";
const juce::String kTipLearnDefault =
    "Loop a noise-only segment in your DAW and click\nto capture a noise "
    "profile.";
const juce::String kTipResetProfile = "Reset noise profile";
const juce::String kTipResetProfileClear =
    "Reset noise profile and clear learned data.";
const juce::String kTipResetProfileUnavailable =
    "There is no profile to reset. Learn a profile or enable Adaptive first.";
const juce::String kTipAdaptiveMethod =
    "Choose how Adaptive estimates the noise floor:\n"
    "Speech-focused (SPP-MMSE): considers speech presence when updating.\n"
    "Slow, steady (Brandt): uses a longer rolling statistical history.\n"
    "Balanced tracking (Martin): follows spectral minima over time.";
// Reduction faders
const juce::String kTipMasterReduction =
    "Adjust noise reduction level in decibels\n(0 to 40 dB across all "
    "bands).";
const juce::String kTipTonalReduction =
    "Adjust reduction level for tonal noise components\n(0 to 40 dB for "
    "harmonic peaks).";
const juce::String kTipDelta =
    "Listen to removed noise signal\n(residual audio output).";
const juce::String kTipBypass =
    "Soft bypass plugin processing\nwith smooth wet/dry fade transition.";
const juce::String kTipLinkUnlinkedDisabled =
    "Unlinking tonal reduction is disabled in Standalone Adaptive "
    "mode\n(requires a captured manual profile to detect tonal peaks).";
const juce::String kTipCurveToggle =
    "Enable per-frequency reduction bias curve.";
const juce::String kTipCurveOverlay =
    "Enable per-frequency reduction bias curve overlay on spectral\n"
    "display. Shapes both broadband and tonal reduction.";
// Threshold faders
const juce::String kTipMasterOffsetUnlinked =
    "Shift broadband noise threshold up or down\nin "
    "decibels (-12 to +12 dB).";
const juce::String kTipTonalOffset =
    "Shift threshold for tonal noise components\nin decibels (-12 to +12 "
    "dB).";
// Advanced panel
const juce::String kTipSmoothing =
    "Apply temporal smoothing across spectral frames\nto reduce musical "
    "noise bubbling artifacts.";
const juce::String kTipSmoothing2D =
    "Adjust patch-based similarity filtering strength\nto "
    "control patch-based (NLM / DFTT / BM3D) smoothing.";
const juce::String kTipMasking =
    "Keep low-level details that loud content masks.\nHigher values "
    "protect more of the residual floor.";
const juce::String kTipWhitening =
    "Spectral whitening factor to equalize\nthe residual noise floor "
    "spectrum.";
const juce::String kTipReductionBlend =
    "Blend bin-wise reduction with neighborhood-wide gains\nto tame "
    "musical noise; high values keep residual noise bursts.";
const juce::String kTipMorphingUnavailable =
    "Profile Morphing is available in Advanced Controls\nwhen a manual "
    "learned noise profile is captured.";

} // namespace

namespace {

// Frame-size entries share one parameter write path.
void addFrameSizeMenuItems(juce::PopupMenu& menu,
                           juce::RangedAudioParameter* frameSizeParam,
                           int frameSizeIdx, bool lowLatencyVal) {
  // Guidance baked into the labels: PopupMenu items have no tooltip
  // support, and latency in ms is sample-rate independent (always 2x).
  static constexpr std::array<const char*, 5> kFrameSizeNames = {
      "23 ms - preserve transients", "32 ms - general purpose",
      "46 ms - stationary noise, default", "64 ms - tonal detail: hum & whine",
      "93 ms - max resolution"};
  for (int i = 0; i < 5; ++i) {
    juce::PopupMenu::Item frameItem(kFrameSizeNames[i]);
    frameItem.itemID = 3 + i;
    frameItem.isTicked = (i == frameSizeIdx);
    frameItem.isEnabled = !lowLatencyVal;
    frameItem.action = [frameSizeParam, i]() {
      if (frameSizeParam != nullptr) {
        frameSizeParam->beginChangeGesture();
        frameSizeParam->setValueNotifyingHost(static_cast<float>(i) / 4.0f);
        frameSizeParam->endChangeGesture();
      }
    };
    menu.addItem(frameItem);
  }
}

// Learn button appearance follows capture / profile state.
void styleLearnButton(juce::TextButton& btn, bool isLearning, bool hasProfile) {
  if (isLearning) {
    btn.setButtonText("Learning...");
    btn.setColour(
        juce::TextButton::buttonColourId,
        juce::Colour(
            NoiseRepellentLookAndFeel::kColorLearnActive)); // Active Learning
                                                            // Red
    btn.setColour(juce::TextButton::buttonOnColourId,
                  juce::Colour(NoiseRepellentLookAndFeel::kColorLearnActive));
    btn.setTooltip(kTipLearnCapturing);
  } else if (hasProfile) {
    btn.setButtonText("+ Learn");
    btn.setColour(juce::TextButton::buttonColourId,
                  NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btn.setColour(juce::TextButton::buttonOnColourId,
                  NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btn.setTooltip(kTipLearnAccumulate);
  } else {
    btn.setButtonText("Learn Noise");
    btn.setColour(
        juce::TextButton::buttonColourId,
        juce::Colour(
            NoiseRepellentLookAndFeel::kColorLearnCTA)); // Prominent Crimson
                                                         // Red CTA
    btn.setColour(
        juce::TextButton::buttonOnColourId,
        juce::Colour(
            NoiseRepellentLookAndFeel::kColorLearnActive)); // Active Learning
                                                            // Red
    btn.setTooltip(kTipLearnDefault);
  }
  btn.repaint();
}

void styleAdaptiveButton(juce::TextButton& btn, bool isAdaptive) {
  if (isAdaptive) {
    btn.setColour(juce::TextButton::buttonColourId,
                  NoiseRepellentLookAndFeel::kColorDenoising);
    btn.setColour(juce::TextButton::buttonOnColourId,
                  NoiseRepellentLookAndFeel::kColorDenoising);
  } else {
    btn.setColour(juce::TextButton::buttonColourId,
                  juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
    btn.setColour(juce::TextButton::buttonOnColourId,
                  NoiseRepellentLookAndFeel::kColorDenoising);
  }
  btn.repaint();
}

// Status label (HUD overlay on FFT spectrum chart)
void setProfileStatusText(juce::Label& lbl, bool isBypassed, bool isDelta,
                          bool isLearning, bool isAdaptive, bool hasProfile) {
  if (isBypassed) {
    lbl.setText("STATUS: BYPASSED", juce::dontSendNotification);
    lbl.setColour(juce::Label::textColourId,
                  juce::Colour(NoiseRepellentLookAndFeel::kColorInactiveText));
  } else if (isDelta) {
    lbl.setText("STATUS: DELTA (RESIDUAL LISTEN)", juce::dontSendNotification);
    lbl.setColour(juce::Label::textColourId,
                  NoiseRepellentLookAndFeel::kColorTonalPeaks);
  } else if (isLearning) {
    lbl.setText("STATUS: CAPTURING NOISE PROFILE...",
                juce::dontSendNotification);
    lbl.setColour(juce::Label::textColourId,
                  juce::Colour(NoiseRepellentLookAndFeel::kColorLearnActive));
  } else if (isAdaptive && hasProfile) {
    lbl.setText("STATUS: HYBRID (ADAPTIVE + PROFILE)",
                juce::dontSendNotification);
    lbl.setColour(juce::Label::textColourId,
                  NoiseRepellentLookAndFeel::kColorNoiseProfile);
  } else if (isAdaptive) {
    lbl.setText("STATUS: CONTINUOUS ADAPTIVE", juce::dontSendNotification);
    lbl.setColour(juce::Label::textColourId,
                  NoiseRepellentLookAndFeel::kColorDenoising);
  } else if (hasProfile) {
    lbl.setText("STATUS: STATIONARY NOISE PROFILE", juce::dontSendNotification);
    lbl.setColour(juce::Label::textColourId,
                  NoiseRepellentLookAndFeel::kColorNoiseProfile);
  } else {
    lbl.setText("STATUS: NO PROFILE (PASS-THROUGH)",
                juce::dontSendNotification);
    lbl.setColour(juce::Label::textColourId,
                  juce::Colour(NoiseRepellentLookAndFeel::kColorInactiveText));
  }
}

} // namespace

// Menu construction without showing it: PopupMenu building needs no display,
// so Linux tests can exercise it (only showMenuAsync requires X11).
juce::PopupMenu::Item makeToggleItem(int itemID, const juce::String& label,
                                     juce::AudioProcessorParameter* param,
                                     bool ticked,
                                     const std::function<void()>& relayout) {
  juce::PopupMenu::Item item(label);
  item.itemID = itemID;
  item.isTicked = ticked;
  item.action = [param, ticked, relayout]() {
    if (param == nullptr)
      return;
    param->beginChangeGesture();
    param->setValueNotifyingHost(ticked ? 0.0f : 1.0f);
    param->endChangeGesture();
    if (relayout != nullptr)
      relayout();
  };
  return item;
}

juce::PopupMenu buildEditorPreferencesMenu(
    NoiseRepellentAudioProcessor& processor, juce::LookAndFeel* menuLookAndFeel,
    const std::function<void()>& relayout) {
  auto* tooltipParam = processor.getAPVTS().getParameter("show_tooltips");
  auto* hudParam = processor.getAPVTS().getParameter("show_hud");
  auto* frameSizeParam = processor.getAPVTS().getParameter("frame_size");
  auto* lowLatencyParam = processor.getAPVTS().getParameter("low_latency");

  bool tooltipVal =
      (tooltipParam != nullptr && tooltipParam->getValue() > 0.5f);
  bool hudVal = (hudParam != nullptr && hudParam->getValue() > 0.5f);
  bool lowLatencyVal =
      (lowLatencyParam != nullptr && lowLatencyParam->getValue() > 0.5f);
  int frameSizeIdx = 2;
  if (const auto* choice =
          dynamic_cast<juce::AudioParameterChoice*>(frameSizeParam))
    frameSizeIdx = choice->getIndex();

  juce::PopupMenu menu;

  menu.addItem(
      makeToggleItem(1, "Show Tooltips", tooltipParam, tooltipVal, relayout));
  menu.addItem(
      makeToggleItem(2, "Show HUD Status", hudParam, hudVal, relayout));
  menu.addSeparator();
  menu.addItem(makeToggleItem(8, "Low Latency (512, live use)", lowLatencyParam,
                              lowLatencyVal, {}));
  menu.addSeparator();
  menu.addSectionHeader("FRAME SIZE (BIGGER = FINER DETAIL + MORE LATENCY)");
  addFrameSizeMenuItems(menu, frameSizeParam, frameSizeIdx, lowLatencyVal);

  if (menuLookAndFeel != nullptr) {
    menu.setLookAndFeel(menuLookAndFeel);
  }
  return menu;
}

// Dialog text assembly is pure string building; only showing the alert
// window needs a display.
juce::PopupMenu buildAdaptiveMethodMenu(int currentMethod,
                                        juce::LookAndFeel* menuLookAndFeel) {
  juce::PopupMenu menu;
  menu.addSectionHeader("ADAPTIVE ESTIMATION METHOD");
  menu.addItem(1, "Speech-Focused (SPP-MMSE)", true, currentMethod == 1);
  menu.addItem(2, "Slow, Steady (Brandt)", true, currentMethod == 2);
  menu.addItem(3, "Balanced Tracking (Martin)", true, currentMethod == 3);
  if (menuLookAndFeel != nullptr) {
    menu.setLookAndFeel(menuLookAndFeel);
  }
  return menu;
}

juce::String buildAboutMessage() {
#ifdef JucePlugin_VersionString
  const juce::String pluginVersion = JucePlugin_VersionString;
#else
  const juce::String pluginVersion = "dev";
#endif
  const juce::String builtAgainst = SPECBLEACH_VERSION_STRING;
  const juce::String expectedBanner =
      juce::String("libspecbleach ") + builtAgainst;
  const char* runtime = specbleach_get_version_string();
  juce::String message;
  message << "Noise Repellent " << pluginVersion << " (GPL-3.0-or-later)\n\n"
          << "Built against libspecbleach " << builtAgainst << "\n"
          << "Loaded libspecbleach: "
          << (runtime != nullptr ? runtime : "unknown"); // LCOV_EXCL_BR_LINE
  if (runtime != nullptr &&
      juce::String(runtime) != expectedBanner) // LCOV_EXCL_BR_LINE
    message << "\n\nNote: the loaded library differs from the build version.";
  return message;
}

NoiseRepellentAudioProcessorEditor::NoiseRepellentAudioProcessorEditor(
    NoiseRepellentAudioProcessor& p)
    : AudioProcessorEditor(&p), audioProcessor(p),
      spectralVisualizer(p) { // LCOV_EXCL_BR_LINE
  setLookAndFeel(&customLookAndFeel);

  // 1. Native Resizable Window Settings
  setResizable(true, true);
  setResizeLimits(840, 480, 1600, 1000);
  setSize(940, 560);

  // Brand Header (accessible button opening the About box with dependency
  // versions)
  header.brandButton.setName(NoiseRepellentLookAndFeel::kBrandButtonName);
  header.brandButton.setButtonText("NOISE REPELLENT");
  header.brandButton.setTooltip(kTipAbout);
  header.brandButton.setMouseCursor(juce::MouseCursor::PointingHandCursor);
  header.brandButton.onClick = [this]() { showAboutBox(); };
  addAndMakeVisible(header.brandButton);

  // Preferences Header Button
  addAndMakeVisible(header.btnPreferences);
  header.btnPreferences.setTooltip(kTipPreferences);
  header.btnPreferences.setColour(juce::TextButton::buttonColourId,
                                  juce::Colour(0xff3a4150));
  header.btnPreferences.onClick = [this]() { showPreferencesMenu(); };

  // Algorithm Mode Label & Dropdown
  header.lblAlgoHeader.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  header.lblAlgoHeader.setColour(juce::Label::textColourId,
                                 NoiseRepellentLookAndFeel::kColorDenoising);
  header.lblAlgoHeader.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(header.lblAlgoHeader);

  header.comboAlgoMode.addItemList(
      {"Temporal (Fast / Low CPU)", "Patch-Based NLM (Texture & Detail)",
       "NLM + DFTT Refinement (Tonal / Difficult)",
       "Collaborative BM3D (Smooth / Steady)"},
      1);
  addAndMakeVisible(header.comboAlgoMode);
  header.comboAlgoMode.onChange = [this]() { updateLayout(); };

  // Header Action Toggles
  isAdvancedVisible =
      audioProcessor.getAPVTS().getRawParameterValue("show_advanced")->load() >
      0.5f;
  header.btnAdvancedToggle.setClickingTogglesState(true);
  addAndMakeVisible(header.btnAdvancedToggle);
  header.btnAdvancedToggle.setColour(
      juce::TextButton::buttonColourId,
      juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  header.btnAdvancedToggle.setColour(
      juce::TextButton::buttonOnColourId,
      NoiseRepellentLookAndFeel::kColorNoiseProfile);
  header.btnAdvancedToggle.onClick = [this]() {
    isAdvancedVisible = header.btnAdvancedToggle.getToggleState();
    updateLayout();
  };

  header.btnDelta.setClickingTogglesState(true);
  header.btnDelta.setColour(
      juce::TextButton::buttonColourId,
      juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  header.btnDelta.setColour(juce::TextButton::buttonOnColourId,
                            NoiseRepellentLookAndFeel::kColorNoiseProfile);
  addAndMakeVisible(header.btnDelta);
  header.btnDelta.onClick = [this]() { updateLayout(); };

  header.btnBypass.setClickingTogglesState(true);
  header.btnBypass.setColour(
      juce::TextButton::buttonColourId,
      juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  header.btnBypass.setColour(juce::TextButton::buttonOnColourId,
                             NoiseRepellentLookAndFeel::kColorNoiseProfile);
  addAndMakeVisible(header.btnBypass);
  header.btnBypass.onClick = [this]() { updateLayout(); };

  // Beginner Module 1: Noise Profile
  addAndMakeVisible(profile.groupProfile);
  profile.groupProfile.setText("NOISE PROFILE");
  profile.groupProfile.setColour(juce::GroupComponent::outlineColourId,
                                 NoiseRepellentLookAndFeel::kColorPanelBorder);
  profile.groupProfile.setColour(juce::GroupComponent::textColourId,
                                 NoiseRepellentLookAndFeel::kColorFineTuning);
  profile.groupProfile.setInterceptsMouseClicks(false, false);

  profile.btnLearn.setClickingTogglesState(true);
  addAndMakeVisible(profile.btnLearn);
  profile.btnLearn.setColour(
      juce::TextButton::buttonColourId,
      juce::Colour(
          NoiseRepellentLookAndFeel::kColorLearnCTA)); // Prominent Studio
                                                       // Crimson Red CTA
  profile.btnLearn.setColour(
      juce::TextButton::buttonOnColourId,
      juce::Colour(
          NoiseRepellentLookAndFeel::kColorLearnActive)); // Active Learning Red
  profile.btnLearn.onClick = [this]() { updateProfileStatus(); };

  profile.btnAdaptiveNoise.setClickingTogglesState(true);
  profile.btnAdaptiveNoise.setButtonText("Adaptive");
  profile.btnAdaptiveNoise.setTooltip(kTipAdaptive);
  addAndMakeVisible(profile.btnAdaptiveNoise);
  profile.btnAdaptiveNoise.setColour(
      juce::TextButton::buttonColourId,
      juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  profile.btnAdaptiveNoise.setColour(
      juce::TextButton::buttonOnColourId,
      NoiseRepellentLookAndFeel::kColorDenoising);
  profile.btnAdaptiveNoise.onClick = [this]() {
    updateLayout();
    updateProfileStatus();
  };

  addAndMakeVisible(profile.btnAdaptiveArrow);
  profile.btnAdaptiveArrow.setTooltip(kTipAdaptiveMethod);
  profile.btnAdaptiveArrow.setColour(juce::TextButton::buttonColourId,
                                     juce::Colour(0xff353b48));
  profile.btnAdaptiveArrow.onClick = [this]() { showAdaptiveMethodMenu(); };

  addAndMakeVisible(profile.btnResetProfile);
  profile.btnResetProfile.setTooltip(kTipResetProfile);
  profile.btnResetProfile.onClick = [this]() { handleResetProfileClick(); };

  // Profile Status Label (HUD Overlay on Chart)
  addAndMakeVisible(profile.lblProfileStatus);
  profile.lblProfileStatus.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  profile.lblProfileStatus.setColour(
      juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
  profile.lblProfileStatus.setJustificationType(
      juce::Justification::centredRight);
  profile.lblProfileStatus.setText("STATUS: STATIONARY NOISE PROFILE",
                                   juce::dontSendNotification);

  // Beginner Module 2: Denoising & Resizable Canvas
  addAndMakeVisible(denoising.groupDenoising);
  denoising.groupDenoising.setText("");
  denoising.groupDenoising.setColour(
      juce::GroupComponent::outlineColourId,
      NoiseRepellentLookAndFeel::kColorPanelBorder);
  denoising.groupDenoising.setColour(
      juce::GroupComponent::textColourId,
      NoiseRepellentLookAndFeel::kColorFineTuning);
  denoising.groupDenoising.setInterceptsMouseClicks(false, false);

  denoising.lblReductionHeader.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  denoising.lblReductionHeader.setColour(
      juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
  denoising.lblReductionHeader.setJustificationType(
      juce::Justification::centred);
  addAndMakeVisible(denoising.lblReductionHeader);

  addAndMakeVisible(denoising.btnLink);
  // Attachment must be created before onClick is assigned: ButtonAttachment
  // calls setToggleState(sendNotification) during construction to sync from
  // APVTS, which would fire onClick and reset reduction params to default.
  attachments.attachLink = std::make_unique<ButtonAttachment>(
      audioProcessor.getAPVTS(), "link_reduction", denoising.btnLink);
  denoising.btnLink.setButtonText(
      denoising.btnLink.getToggleState() ? "Linked" : "Unlinked");
  denoising.btnLink.onClick = [this]() { handleLinkToggled(); };

  denoising.btnCurveToggle.setClickingTogglesState(true);
  denoising.btnCurveToggle.setTooltip(kTipCurveToggle);
  addAndMakeVisible(denoising.btnCurveToggle);
  denoising.btnCurveToggle.onClick = [this]() { updateLayout(); };

  denoising.btnResetCurve.setTooltip(kTipResetCurve);
  addAndMakeVisible(denoising.btnResetCurve);
  denoising.btnResetCurve.onClick = [this]() {
    audioProcessor.curve.resetNodes();
    spectralVisualizer.repaint();
  };

  denoising.lblMasterRed.setText("BROADBAND", juce::dontSendNotification);
  denoising.lblMasterRed.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  denoising.lblMasterRed.setColour(juce::Label::textColourId,
                                   NoiseRepellentLookAndFeel::kColorDenoising);
  denoising.lblMasterRed.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(denoising.lblMasterRed);

  denoising.lblTonalRed.setText("TONAL", juce::dontSendNotification);
  denoising.lblTonalRed.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  denoising.lblTonalRed.setColour(juce::Label::textColourId,
                                  NoiseRepellentLookAndFeel::kColorDenoising);
  denoising.lblTonalRed.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(denoising.lblTonalRed);

  denoising.sliderMasterRed.setColour(
      juce::Slider::rotarySliderFillColourId,
      NoiseRepellentLookAndFeel::kColorDenoising);
  denoising.sliderTonalRed.setColour(
      juce::Slider::rotarySliderFillColourId,
      NoiseRepellentLookAndFeel::kColorDenoising);

  addAndMakeVisible(denoising.sliderMasterRed);
  addAndMakeVisible(denoising.sliderTonalRed);

  addAndMakeVisible(spectralVisualizer);

  // Collapsible Module 3: Advanced Controls
  addAndMakeVisible(advanced.groupAdvanced);
  advanced.groupAdvanced.setText("ADVANCED CONTROLS");
  advanced.groupAdvanced.setColour(
      juce::GroupComponent::outlineColourId,
      NoiseRepellentLookAndFeel::kColorPanelBorder);
  advanced.groupAdvanced.setColour(juce::GroupComponent::textColourId,
                                   NoiseRepellentLookAndFeel::kColorFineTuning);
  advanced.groupAdvanced.setInterceptsMouseClicks(false, false);

  advanced.comboMethod.addItemList(
      {"Speech-Focused (SPP-MMSE)", "Slow, Steady (Brandt)",
       "Balanced Tracking (Martin)"},
      1);
  addAndMakeVisible(advanced.comboMethod);

  advanced.lblMethod.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  advanced.lblMethod.setColour(juce::Label::textColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.lblMethod.setJustificationType(juce::Justification::left);
  addAndMakeVisible(advanced.lblMethod);

  profile.lblOffset.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  profile.lblOffset.setColour(juce::Label::textColourId,
                              NoiseRepellentLookAndFeel::kColorDenoising);
  profile.lblOffset.setJustificationType(juce::Justification::centred);
  profile.sliderOffset.setColour(juce::Slider::rotarySliderFillColourId,
                                 NoiseRepellentLookAndFeel::kColorNoiseProfile);
  addAndMakeVisible(profile.sliderOffset);
  addAndMakeVisible(profile.lblOffset);

  profile.lblMasterOffset.setText("BROADBAND", juce::dontSendNotification);
  profile.lblMasterOffset.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  profile.lblMasterOffset.setColour(
      juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
  profile.lblMasterOffset.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(profile.lblMasterOffset);

  profile.lblTonalOffset.setText("TONAL", juce::dontSendNotification);
  profile.lblTonalOffset.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  profile.lblTonalOffset.setColour(juce::Label::textColourId,
                                   NoiseRepellentLookAndFeel::kColorTonalPeaks);
  profile.lblTonalOffset.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(profile.lblTonalOffset);

  profile.sliderTonalOffset.setColour(
      juce::Slider::rotarySliderFillColourId,
      NoiseRepellentLookAndFeel::kColorTonalPeaks);

  addAndMakeVisible(profile.sliderTonalOffset);

  addAndMakeVisible(profile.btnLinkOffset);
  // Attachment must be created before onClick is assigned (same reasoning as
  // the reduction link button above).
  attachments.attachLinkOffset = std::make_unique<ButtonAttachment>(
      audioProcessor.getAPVTS(), "link_threshold_offset",
      profile.btnLinkOffset);
  profile.btnLinkOffset.setButtonText(
      profile.btnLinkOffset.getToggleState() ? "Linked" : "Unlinked");
  profile.btnLinkOffset.onClick = [this]() { handleLinkOffsetToggled(); };

  advanced.lblAggressiveness.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  advanced.lblAggressiveness.setColour(
      juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.lblAggressiveness.setJustificationType(juce::Justification::centred);
  advanced.sliderAggressiveness.setColour(
      juce::Slider::rotarySliderFillColourId,
      NoiseRepellentLookAndFeel::kColorDenoising);
  addAndMakeVisible(advanced.sliderAggressiveness);
  addAndMakeVisible(advanced.lblAggressiveness);

  advanced.sliderSmoothing.setColour(
      juce::Slider::rotarySliderFillColourId,
      NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.sliderMasking.setColour(juce::Slider::rotarySliderFillColourId,
                                   NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.sliderWhitening.setColour(
      juce::Slider::rotarySliderFillColourId,
      NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.sliderReductionBlend.setColour(
      juce::Slider::rotarySliderFillColourId,
      NoiseRepellentLookAndFeel::kColorDenoising);

  addAndMakeVisible(advanced.sliderSmoothing);
  addAndMakeVisible(advanced.sliderMasking);
  addAndMakeVisible(advanced.sliderWhitening);
  addAndMakeVisible(advanced.sliderReductionBlend);

  advanced.lblSmoothing.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  advanced.lblSmoothing.setColour(juce::Label::textColourId,
                                  NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.lblSmoothing.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(advanced.lblSmoothing);

  advanced.lblMasking.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  advanced.lblMasking.setColour(juce::Label::textColourId,
                                NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.lblMasking.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(advanced.lblMasking);

  advanced.lblWhitening.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  advanced.lblWhitening.setColour(juce::Label::textColourId,
                                  NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.lblWhitening.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(advanced.lblWhitening);

  advanced.lblReductionBlend.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  advanced.lblReductionBlend.setColour(
      juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
  advanced.lblReductionBlend.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(advanced.lblReductionBlend);

  // Advanced Controls Toggle Button
  header.btnAdvancedToggle.setClickingTogglesState(true);
  header.btnAdvancedToggle.setButtonText("ADVANCED");
  header.btnAdvancedToggle.setTooltip(kTipAdvancedToggle);
  addAndMakeVisible(header.btnAdvancedToggle);
  header.btnAdvancedToggle.onClick = [this]() {
    isAdvancedVisible = header.btnAdvancedToggle.getToggleState();
    updateLayout();
  };

  // Parameter APVTS Attachments
  auto& apvts = audioProcessor.getAPVTS();

  attachments.attachAlgoMode = std::make_unique<ComboBoxAttachment>(
      apvts, "algorithm_mode", header.comboAlgoMode);
  attachments.attachLearn = std::make_unique<ButtonAttachment>(
      apvts, "learn_noise", profile.btnLearn);
  attachments.attachAdaptive = std::make_unique<ButtonAttachment>(
      apvts, "adaptive_noise", profile.btnAdaptiveNoise);
  attachments.attachCurveToggle = std::make_unique<ButtonAttachment>(
      apvts, "reduction_curve_enabled", denoising.btnCurveToggle);
  attachments.attachDelta = std::make_unique<ButtonAttachment>(
      apvts, "residual_listen", header.btnDelta);
  attachments.attachBypass =
      std::make_unique<ButtonAttachment>(apvts, "bypass", header.btnBypass);
  attachments.attachShowAdvanced = std::make_unique<ButtonAttachment>(
      apvts, "show_advanced", header.btnAdvancedToggle);

  attachments.attachMethod = std::make_unique<ComboBoxAttachment>(
      apvts, "adaptive_method", advanced.comboMethod);

  attachments.attachMasterRed = std::make_unique<SliderAttachment>(
      apvts, "reduction_amount", denoising.sliderMasterRed);
  attachments.attachTonalRed = std::make_unique<SliderAttachment>(
      apvts, "tonal_reduction", denoising.sliderTonalRed);
  attachments.attachOffset = std::make_unique<SliderAttachment>(
      apvts, "noise_profile_offset", profile.sliderOffset);
  attachments.attachTonalOffset = std::make_unique<SliderAttachment>(
      apvts, "tonal_noise_profile_offset", profile.sliderTonalOffset);
  attachments.attachSmoothing = std::make_unique<SliderAttachment>(
      apvts, "smoothing_factor", advanced.sliderSmoothing);
  attachments.attachMasking = std::make_unique<SliderAttachment>(
      apvts, "masking_depth", advanced.sliderMasking);
  attachments.attachWhitening = std::make_unique<SliderAttachment>(
      apvts, "whitening_factor", advanced.sliderWhitening);
  attachments.attachAggressiveness = std::make_unique<SliderAttachment>(
      apvts, "aggressiveness", advanced.sliderAggressiveness);
  attachments.attachReductionBlend = std::make_unique<SliderAttachment>(
      apvts, "reduction_blend", advanced.sliderReductionBlend);

  // Control Tooltip Descriptions
  header.btnPreferences.setTooltip(kTipPreferencesMenu);
  header.comboAlgoMode.setTooltip(kTipAlgoMode);
  header.lblAlgoHeader.setTooltip(kTipAlgoMode);
  header.btnAdvancedToggle.setTooltip(kTipAdvancedShow);
  profile.btnLearn.setTooltip(kTipLearnInitial);
  profile.btnResetProfile.setTooltip(kTipResetProfileClear);
  profile.btnAdaptiveNoise.setTooltip(kTipAdaptive);
  profile.btnAdaptiveArrow.setTooltip(kTipAdaptiveMethod);
  advanced.comboMethod.setTooltip(kTipAdaptiveMethod);
  advanced.lblMethod.setTooltip(kTipAdaptiveMethod);

  denoising.sliderMasterRed.setTooltip(kTipMasterReduction);
  denoising.lblMasterRed.setTooltip(kTipMasterReduction);
  denoising.lblReductionHeader.setTooltip(kTipMasterReduction);
  denoising.sliderTonalRed.setTooltip(kTipTonalReduction);
  denoising.lblTonalRed.setTooltip(kTipTonalReduction);

  profile.sliderOffset.setTooltip(kTipThreshold);
  profile.lblOffset.setTooltip(kTipThreshold);
  advanced.sliderAggressiveness.setTooltip(kTipAggressiveness);
  advanced.lblAggressiveness.setTooltip(kTipAggressiveness);
  header.btnDelta.setTooltip(kTipDelta);
  header.btnBypass.setTooltip(kTipBypass);

  advanced.sliderSmoothing.setTooltip(kTipSmoothing);
  advanced.sliderMasking.setTooltip(kTipMasking);
  advanced.sliderWhitening.setTooltip(kTipWhitening);
  advanced.sliderReductionBlend.setTooltip(kTipReductionBlend);
  advanced.lblReductionBlend.setTooltip(kTipReductionBlend);
  denoising.btnLink.setTooltip(kTipLink);
  denoising.btnCurveToggle.setTooltip(kTipCurveOverlay);
  denoising.btnResetCurve.setTooltip(kTipResetCurve);

  // Initial Layout update
  isAdvancedVisible = header.btnAdvancedToggle.getToggleState();
  updateLayout();
  updateSliderLabels();
  // Seed the offline-activity baseline so a stale kOffline flag (e.g. from
  // Audacity's load-time setupProcessing) doesn't flash the overlay on open.
  lastSeenOfflineBlockCount = audioProcessor.status.getNonRealtimeBlockCount();
  startTimerHz(30);

  // Initial default instruction text
  footerTooltipLabel.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeTooltip, juce::Font::plain));
  footerTooltipLabel.setColour(
      juce::Label::textColourId,
      juce::Colour(NoiseRepellentLookAndFeel::kColorFooterText));
  footerTooltipLabel.setJustificationType(juce::Justification::left);
  addAndMakeVisible(footerTooltipLabel);

  if (const auto* showTooltipsParam =
          audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
      showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f) {
    footerTooltipLabel.setText(kTipDefaultInstruction,
                               juce::dontSendNotification);
  }

  // Register parameter listener for show_tooltips
  audioProcessor.getAPVTS().addParameterListener("show_tooltips", this);

  // Register mouse listener AFTER all child components are added
  addMouseListener(this, true);
}

void NoiseRepellentAudioProcessorEditor::
    showPreferencesMenu() { // LCOV_EXCL_START async menus need a real display
  buildEditorPreferencesMenu(audioProcessor, &getLookAndFeel(),
                             [this]() { updateLayout(); })
      .showMenuAsync(juce::PopupMenu::Options().withTargetComponent(
          &header.btnPreferences));
} // LCOV_EXCL_STOP

void NoiseRepellentAudioProcessorEditor::
    showAdaptiveMethodMenu() { // LCOV_EXCL_START async menus need a real
                               // display
  buildAdaptiveMethodMenu(advanced.comboMethod.getSelectedId(),
                          &getLookAndFeel())
      .showMenuAsync(juce::PopupMenu::Options().withTargetComponent(
                         &profile.btnAdaptiveArrow),
                     [this](int result) { applyAdaptiveMethodChoice(result); });
} // LCOV_EXCL_STOP

void NoiseRepellentAudioProcessorEditor::applyAdaptiveMethodChoice(int result) {
  if (result >= 1 && result <= 3) {
    int methodIdx = result - 1;
    advanced.comboMethod.setSelectedId(result, juce::sendNotification);
    if (auto* p = audioProcessor.getAPVTS().getParameter("adaptive_method"))
      p->setValueNotifyingHost(static_cast<float>(methodIdx) / 2.0f);
    if (auto* p = audioProcessor.getAPVTS().getParameter("adaptive_noise"))
      p->setValueNotifyingHost(1.0f);
    updateLayout();
    updateProfileStatus();
  }
}

void NoiseRepellentAudioProcessorEditor::handleResetProfileClick() {
  audioProcessor.resetNoiseProfile();
  if (auto* pLearn = audioProcessor.getAPVTS().getParameter("learn_noise"))
    pLearn->setValueNotifyingHost(0.0f);
  if (auto* pAdapt = audioProcessor.getAPVTS().getParameter("adaptive_noise"))
    pAdapt->setValueNotifyingHost(0.0f);
  updateLayout();
  updateProfileStatus();
}

void NoiseRepellentAudioProcessorEditor::handleLinkToggled() {
  bool isLinked = denoising.btnLink.getToggleState();
  denoising.btnLink.setButtonText(isLinked ? "Linked" : "Unlinked");
  if (isLinked) {
    // Reset reduction parameters to default when relinked
    if (auto* pMaster =
            audioProcessor.getAPVTS().getParameter("reduction_amount"))
      pMaster->setValueNotifyingHost(pMaster->getDefaultValue());
    if (auto* pTonal =
            audioProcessor.getAPVTS().getParameter("tonal_reduction"))
      pTonal->setValueNotifyingHost(pTonal->getDefaultValue());
  }
  updateLayout();
}

void NoiseRepellentAudioProcessorEditor::handleLinkOffsetToggled() {
  bool isLinked = profile.btnLinkOffset.getToggleState();
  profile.btnLinkOffset.setButtonText(isLinked ? "Linked" : "Unlinked");
  if (isLinked) {
    // Reset threshold offset parameters to default when relinked
    if (auto* pMaster =
            audioProcessor.getAPVTS().getParameter("noise_profile_offset"))
      pMaster->setValueNotifyingHost(pMaster->getDefaultValue());
    if (auto* pTonal = audioProcessor.getAPVTS().getParameter(
            "tonal_noise_profile_offset"))
      pTonal->setValueNotifyingHost(pTonal->getDefaultValue());
  }
  updateLayout();
}

void NoiseRepellentAudioProcessorEditor::mouseEnter(
    const juce::MouseEvent& event) {
  if (const auto* showTooltipsParam =
          audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
      showTooltipsParam != nullptr && showTooltipsParam->load() < 0.5f) {
    footerTooltipLabel.setText({}, juce::dontSendNotification);
    return;
  }

  juce::Component* comp = event.originalComponent;
  if (comp == nullptr || comp == this || comp == &profile.groupProfile ||
      comp == &denoising.groupDenoising || comp == &advanced.groupAdvanced) {
    comp = getComponentAt(event.getEventRelativeTo(this).getPosition());
  }

  while (comp != nullptr && comp != this) {
    if (auto* client = dynamic_cast<juce::TooltipClient*>(comp)) {
      auto tip = client->getTooltip();
      if (tip.isNotEmpty()) {
        footerTooltipLabel.setText(tip, juce::dontSendNotification);
        return;
      }
    }
    comp = comp->getParentComponent();
  }

  // Default instruction when hovering over background / non-tooltip components
  footerTooltipLabel.setText(kTipDefaultInstruction,
                             juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::showAboutBox()
    const { // LCOV_EXCL_START modal needs a real display
  juce::AlertWindow::showMessageBoxAsync(juce::MessageBoxIconType::InfoIcon,
                                         "About Noise Repellent",
                                         buildAboutMessage(), "OK");
} // LCOV_EXCL_STOP

void NoiseRepellentAudioProcessorEditor::mouseMove(
    const juce::MouseEvent& event) {
  mouseEnter(event);
}

void NoiseRepellentAudioProcessorEditor::mouseExit(
    const juce::MouseEvent& /*event*/) {
  if (const auto* showTooltipsParam =
          audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
      showTooltipsParam != nullptr && showTooltipsParam->load() > 0.5f) {
    footerTooltipLabel.setText(kTipDefaultInstruction,
                               juce::dontSendNotification);
  } else {
    footerTooltipLabel.setText({}, juce::dontSendNotification);
  }
}

NoiseRepellentAudioProcessorEditor::~NoiseRepellentAudioProcessorEditor() {
  audioProcessor.getAPVTS().removeParameterListener("show_tooltips", this);
  cancelPendingUpdate();
  stopTimer();
  setLookAndFeel(nullptr);
}

void NoiseRepellentAudioProcessorEditor::parameterChanged(
    const juce::String& parameterID, float /*newValue*/) {
  if (parameterID == "show_tooltips") {
    tooltipStateDirty.store(true);
  }
}

void NoiseRepellentAudioProcessorEditor::handleAsyncUpdate() {
  if (const auto* showTooltipsParam =
          audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
      showTooltipsParam != nullptr && showTooltipsParam->load() > 0.5f) {
    footerTooltipLabel.setText(kTipDefaultInstruction,
                               juce::dontSendNotification);
  } else {
    footerTooltipLabel.setText({}, juce::dontSendNotification);
  }
  updateLayout();
}

void NoiseRepellentAudioProcessorEditor::updateLayout() {
  bool isLinked = denoising.btnLink.getToggleState();
  bool isCurveOn = denoising.btnCurveToggle.getToggleState();

  denoising.sliderMasterRed.setVisible(true);
  denoising.sliderTonalRed.setVisible(isAdvancedVisible && !isLinked);

  denoising.lblMasterRed.setVisible(isAdvancedVisible && !isLinked);
  denoising.lblTonalRed.setVisible(isAdvancedVisible && !isLinked);

  denoising.btnLink.setVisible(isAdvancedVisible);
  denoising.btnCurveToggle.setVisible(isAdvancedVisible);
  denoising.btnResetCurve.setVisible(isAdvancedVisible && isCurveOn);

  profile.sliderOffset.setVisible(true);
  profile.lblOffset.setVisible(true);

  bool isOffsetLinked =
      !isAdvancedVisible || profile.btnLinkOffset.getToggleState();
  profile.sliderTonalOffset.setVisible(isAdvancedVisible && !isOffsetLinked);
  profile.lblMasterOffset.setVisible(isAdvancedVisible && !isOffsetLinked);
  profile.lblTonalOffset.setVisible(isAdvancedVisible && !isOffsetLinked);
  profile.btnLinkOffset.setVisible(isAdvancedVisible);

  advanced.sliderAggressiveness.setVisible(isAdvancedVisible);
  advanced.lblAggressiveness.setVisible(isAdvancedVisible);

  advanced.groupAdvanced.setVisible(isAdvancedVisible);
  profile.btnLearn.setVisible(true);
  profile.btnAdaptiveNoise.setVisible(true);
  profile.btnAdaptiveArrow.setVisible(isAdvancedVisible);
  profile.btnResetProfile.setVisible(true);
  advanced.comboMethod.setVisible(false);
  advanced.lblMethod.setVisible(false);
  advanced.sliderSmoothing.setVisible(isAdvancedVisible);
  advanced.lblSmoothing.setVisible(isAdvancedVisible);
  advanced.sliderMasking.setVisible(isAdvancedVisible);
  advanced.lblMasking.setVisible(isAdvancedVisible);
  advanced.sliderWhitening.setVisible(isAdvancedVisible);
  advanced.lblWhitening.setVisible(isAdvancedVisible);
  advanced.sliderReductionBlend.setVisible(isAdvancedVisible);
  advanced.lblReductionBlend.setVisible(isAdvancedVisible);

  spectralVisualizer.setAdvancedControlsVisible(isAdvancedVisible);

  header.btnAdvancedToggle.setButtonText(juce::CharPointer_UTF8(
      isAdvancedVisible ? "ADVANCED \xe2\x96\xb2" : "ADVANCED \xe2\x96\xbc"));

  const auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  const auto* showHudParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_hud");

  bool showTooltips =
      (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f);
  bool showHud = (showHudParam == nullptr || showHudParam->load() > 0.5f);

  footerTooltipLabel.setVisible(showTooltips);
  profile.lblProfileStatus.setVisible(showHud);

  resized();
}

void NoiseRepellentAudioProcessorEditor::updateSliderLabels() {
  bool isCompact = getWidth() < 900;
  profile.lblOffset.setText(isCompact ? "THRESH" : "THRESHOLD",
                            juce::dontSendNotification);
  advanced.lblAggressiveness.setText(
      "AGGRESSIVENESS: " +
          juce::String(advanced.sliderAggressiveness.getValue(), 1),
      juce::dontSendNotification);
  advanced.lblSmoothing.setText(
      "SMOOTHING: " +
          juce::String(static_cast<int>(advanced.sliderSmoothing.getValue())) +
          "%",
      juce::dontSendNotification);
  advanced.lblMasking.setText(
      isCompact ? "MASKING: " + // LCOV_EXCL_BR_LINE
                      juce::String(
                          static_cast<int>(advanced.sliderMasking.getValue())) +
                      "%" // LCOV_EXCL_BR_LINE
                : "MASKING PROTECT: " +
                      juce::String(
                          static_cast<int>(advanced.sliderMasking.getValue())) +
                      "%", // LCOV_EXCL_BR_LINE
      juce::dontSendNotification);
  advanced.lblWhitening.setText(
      "WHITENING: " +
          juce::String(static_cast<int>(advanced.sliderWhitening.getValue())) +
          "%",
      juce::dontSendNotification);
  advanced.lblReductionBlend.setText(
      isCompact ? "BLEND: " +
                      juce::String(static_cast<int>(
                          advanced.sliderReductionBlend.getValue())) +
                      "%"
                : "REDUCTION BLEND: " +
                      juce::String(static_cast<int>(
                          advanced.sliderReductionBlend.getValue())) +
                      "%",
      juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::updateProfileStatus() {
  const auto* adaptParam =
      audioProcessor.getAPVTS().getRawParameterValue("adaptive_noise");
  bool isAdaptive = adaptParam != nullptr && adaptParam->load() > 0.5f;
  bool isLearning = profile.btnLearn.getToggleState();
  bool hasProfile = audioProcessor.hasNoiseProfile();

  styleLearnButton(profile.btnLearn, isLearning, hasProfile);
  styleAdaptiveButton(profile.btnAdaptiveNoise, isAdaptive);

  // When bypassed, disable everything except the Bypass button itself
  bool isBypassed = header.btnBypass.getToggleState();
  bool pluginActive = !isBypassed;
  const bool processingAvailable = hasProfile || isAdaptive;
  const bool controlsActive = pluginActive && processingAvailable;

  const bool lowLatency =
      audioProcessor.isLowLatency() ||
      (audioProcessor.getAPVTS().getRawParameterValue("low_latency") !=
           nullptr &&
       audioProcessor.getAPVTS().getRawParameterValue("low_latency")->load() >
           0.5f);

  applyHeaderStates(controlsActive, pluginActive, processingAvailable,
                    lowLatency);
  applyProfileStates(pluginActive, controlsActive, processingAvailable,
                     isAdaptive, hasProfile);
  applyReductionStates(controlsActive, processingAvailable, lowLatency,
                       isAdaptive, hasProfile);
  applyOffsetStates(controlsActive, processingAvailable, lowLatency, isAdaptive,
                    hasProfile);
  applyAdvancedStates(controlsActive, processingAvailable, lowLatency);
  setProfileStatusText(profile.lblProfileStatus, isBypassed,
                       header.btnDelta.getToggleState(), isLearning, isAdaptive,
                       hasProfile);
}

void NoiseRepellentAudioProcessorEditor::applyHeaderStates(
    bool controlsActive, bool pluginActive, bool processingAvailable,
    bool lowLatency) {
  // Header controls
  header.btnDelta.setEnabled(controlsActive);
  header.btnAdvancedToggle.setEnabled(pluginActive);
  header.comboAlgoMode.setEnabled(controlsActive && !lowLatency);
  if (!processingAvailable)
    header.comboAlgoMode.setTooltip(kTipAwaitingProfile);
  else if (lowLatency)
    header.comboAlgoMode.setTooltip(
        "Locked to Temporal in Low Latency mode (causal 1D-only)");
  else
    header.comboAlgoMode.setTooltip(kTipAlgoMode);
  header.lblAlgoHeader.setEnabled(controlsActive);
}

void NoiseRepellentAudioProcessorEditor::applyProfileStates(
    bool pluginActive, bool controlsActive, bool processingAvailable,
    bool isAdaptive, bool hasProfile) {
  // Noise Profile box: enabled whenever plugin is active
  bool profileEnabled = pluginActive;
  profile.groupProfile.setEnabled(pluginActive);
  denoising.groupDenoising.setEnabled(controlsActive);
  profile.btnLearn.setEnabled(profileEnabled);
  profile.btnAdaptiveNoise.setEnabled(profileEnabled);
  profile.btnAdaptiveArrow.setEnabled(profileEnabled);
  profile.btnResetProfile.setEnabled(profileEnabled &&
                                     (hasProfile || isAdaptive));
  profile.btnResetProfile.setTooltip((hasProfile || isAdaptive)
                                         ? kTipResetProfileClear
                                         : kTipResetProfileUnavailable);
  profile.sliderOffset.setEnabled(controlsActive);
  profile.lblOffset.setEnabled(controlsActive);
  if (!processingAvailable) {
    header.btnDelta.setTooltip(kTipAwaitingProfile);
    profile.sliderOffset.setTooltip(kTipAwaitingProfile);
    profile.lblOffset.setTooltip(kTipAwaitingProfile);
  } else {
    header.btnDelta.setTooltip(kTipDelta);
  }

  // Aggressiveness (Profile Morphing) enabled in Noise Profile box when a
  // manual profile exists and Advanced Controls is ON
  bool aggressivenessEnabled = pluginActive && isAdvancedVisible && hasProfile;
  advanced.sliderAggressiveness.setEnabled(aggressivenessEnabled);
  advanced.lblAggressiveness.setEnabled(aggressivenessEnabled);

  if (aggressivenessEnabled) {
    advanced.sliderAggressiveness.setTooltip(kTipAggressiveness);
    advanced.lblAggressiveness.setTooltip(kTipAggressiveness);
  } else {
    advanced.sliderAggressiveness.setTooltip(kTipMorphingUnavailable);
    advanced.lblAggressiveness.setTooltip(kTipMorphingUnavailable);
  }
}

void NoiseRepellentAudioProcessorEditor::applyReductionStates(
    bool controlsActive, bool processingAvailable, bool lowLatency,
    bool isAdaptive, bool hasProfile) {
  // Reduction controls (low-latency 512-sample frame is too coarse for
  // independent tonal/broadband control: unlinking is unavailable)
  bool canUnlink = (!isAdaptive || hasProfile) && !lowLatency;
  bool allowUnlink = controlsActive && canUnlink;
  denoising.btnLink.setEnabled(allowUnlink);

  if (!canUnlink && !denoising.btnLink.getToggleState()) {
    denoising.btnLink.setToggleState(true, juce::dontSendNotification);
    denoising.btnLink.setButtonText("Linked");
    if (auto* linkParam =
            audioProcessor.getAPVTS().getParameter("link_reduction"))
      linkParam->setValueNotifyingHost(1.0f);
  }

  if (!processingAvailable) {
    denoising.btnLink.setTooltip(kTipAwaitingProfile);
    profile.btnLinkOffset.setTooltip(kTipAwaitingProfile);
  } else if (canUnlink) {
    denoising.btnLink.setTooltip(kTipLink);
    profile.btnLinkOffset.setTooltip(kTipLink);
  } else {
    denoising.btnLink.setTooltip(kTipLinkUnlinkedDisabled);
    profile.btnLinkOffset.setTooltip(kTipLinkUnlinkedDisabled);
  }

  denoising.lblReductionHeader.setEnabled(controlsActive);
  denoising.sliderMasterRed.setEnabled(controlsActive);
  denoising.lblMasterRed.setEnabled(controlsActive);

  bool isLinked = denoising.btnLink.getToggleState();
  bool tonalEnabled = controlsActive && !isLinked && allowUnlink;
  denoising.sliderTonalRed.setEnabled(tonalEnabled);
  denoising.lblTonalRed.setEnabled(tonalEnabled);

  const juce::String masterReductionTip =
      isLinked
          ? "Adjust noise reduction level in decibels\n(0 to 40 dB across all "
            "bands)."
          : "Adjust broadband noise reduction level in decibels\n(0 to 40 dB).";
  const juce::String tonalReductionTip =
      "Adjust reduction level for tonal noise components\n(0 to 40 dB).";

  denoising.lblReductionHeader.setTooltip(
      processingAvailable ? masterReductionTip : kTipAwaitingProfile);
  denoising.sliderMasterRed.setTooltip(
      processingAvailable ? masterReductionTip : kTipAwaitingProfile);
  denoising.lblMasterRed.setTooltip(processingAvailable ? masterReductionTip
                                                        : kTipAwaitingProfile);
  denoising.sliderTonalRed.setTooltip(
      processingAvailable ? tonalReductionTip : kTipAwaitingProfile);
  denoising.lblTonalRed.setTooltip(processingAvailable ? tonalReductionTip
                                                       : kTipAwaitingProfile);
}

void NoiseRepellentAudioProcessorEditor::applyOffsetStates(
    bool controlsActive, bool processingAvailable, bool lowLatency,
    bool isAdaptive, bool hasProfile) {
  // Threshold Offset controls
  bool canUnlink = (!isAdaptive || hasProfile) && !lowLatency;
  bool allowUnlink = controlsActive && canUnlink;
  profile.btnLinkOffset.setEnabled(allowUnlink);

  if (!canUnlink && !profile.btnLinkOffset.getToggleState()) {
    profile.btnLinkOffset.setToggleState(true, juce::dontSendNotification);
    profile.btnLinkOffset.setButtonText("Linked");
    if (auto* linkOffsetParam =
            audioProcessor.getAPVTS().getParameter("link_threshold_offset"))
      linkOffsetParam->setValueNotifyingHost(1.0f);
  }

  bool isOffsetLinked = profile.btnLinkOffset.getToggleState();
  bool tonalOffsetEnabled = controlsActive && !isOffsetLinked && allowUnlink;
  profile.sliderTonalOffset.setEnabled(tonalOffsetEnabled);
  profile.lblTonalOffset.setEnabled(tonalOffsetEnabled);

  const juce::String offsetTip =
      isOffsetLinked ? kTipThreshold : kTipMasterOffsetUnlinked;
  profile.lblOffset.setTooltip(processingAvailable ? offsetTip
                                                   : kTipAwaitingProfile);
  profile.sliderOffset.setTooltip(processingAvailable ? offsetTip
                                                      : kTipAwaitingProfile);
  profile.lblMasterOffset.setTooltip(processingAvailable ? offsetTip
                                                         : kTipAwaitingProfile);
  profile.sliderTonalOffset.setTooltip(
      processingAvailable ? kTipTonalOffset : kTipAwaitingProfile);
  profile.lblTonalOffset.setTooltip(processingAvailable ? kTipTonalOffset
                                                        : kTipAwaitingProfile);

  denoising.btnCurveToggle.setEnabled(controlsActive);
  denoising.btnResetCurve.setEnabled(controlsActive &&
                                     denoising.btnCurveToggle.getToggleState());
  denoising.btnCurveToggle.setTooltip(
      processingAvailable ? kTipCurveToggle : kTipAwaitingProfile);
  denoising.btnResetCurve.setTooltip(processingAvailable ? kTipResetCurve
                                                         : kTipAwaitingProfile);
}

void NoiseRepellentAudioProcessorEditor::applyAdvancedStates(
    bool controlsActive, bool processingAvailable, bool lowLatency) {
  // Advanced controls
  advanced.sliderSmoothing.setEnabled(controlsActive);
  advanced.sliderWhitening.setEnabled(controlsActive);
  advanced.sliderReductionBlend.setEnabled(controlsActive);
  advanced.lblSmoothing.setEnabled(controlsActive);
  advanced.lblWhitening.setEnabled(controlsActive);
  advanced.lblReductionBlend.setEnabled(controlsActive);
  advanced.sliderMasking.setEnabled(controlsActive);
  advanced.lblMasking.setEnabled(controlsActive);
  advanced.comboMethod.setEnabled(controlsActive);
  advanced.lblMethod.setEnabled(controlsActive);
  advanced.groupAdvanced.setEnabled(controlsActive);

  // In low latency the advanced panel swaps the Masking slider for the
  // Reduction Blend slider; neither transition changes their values.
  advanced.sliderMasking.setVisible(isAdvancedVisible && !lowLatency);
  advanced.lblMasking.setVisible(isAdvancedVisible && !lowLatency);
  advanced.sliderReductionBlend.setVisible(isAdvancedVisible);
  advanced.lblReductionBlend.setVisible(isAdvancedVisible);
  if (lowLatency != lastAdvancedLowLatency) {
    lastAdvancedLowLatency = lowLatency;
    resized(); // panel column count depends on the active slider pair
  }
  spectralVisualizer.setProcessingAvailable(processingAvailable);

  bool is2D = (header.comboAlgoMode.getSelectedItemIndex() >= 1);
  const juce::String smoothingTip = is2D ? kTipSmoothing2D : kTipSmoothing;

  advanced.sliderSmoothing.setTooltip(
      processingAvailable ? smoothingTip : kTipAwaitingProfile);
  advanced.lblSmoothing.setTooltip(processingAvailable ? smoothingTip
                                                       : kTipAwaitingProfile);
  advanced.sliderMasking.setTooltip(processingAvailable ? kTipMasking
                                                        : kTipAwaitingProfile);
  advanced.lblMasking.setTooltip(processingAvailable ? kTipMasking
                                                     : kTipAwaitingProfile);
  advanced.sliderWhitening.setTooltip(
      processingAvailable ? kTipWhitening : kTipAwaitingProfile);
  advanced.lblWhitening.setTooltip(processingAvailable ? kTipWhitening
                                                       : kTipAwaitingProfile);
}

void NoiseRepellentAudioProcessorEditor::timerCallback() {
  if (tooltipStateDirty.exchange(false)) {
    triggerAsyncUpdate();
  }

  // Audacity pre-initializes VST3 with kOffline at load and only switches
  // to kRealtime when its realtime stack runs, so the raw non-realtime flag
  // is stuck on while the editor sits idle. Treat offline as active only
  // while offline blocks are actually flowing through the processor.
  const uint64_t offlineBlocks =
      audioProcessor.status.getNonRealtimeBlockCount();
  const bool isOffline = audioProcessor.isNonRealtime() &&
                         (offlineBlocks != lastSeenOfflineBlockCount);
  lastSeenOfflineBlockCount = offlineBlocks;
  if (isOffline != wasOfflineRendering) {
    wasOfflineRendering = isOffline;
    repaint();
  }

  updateSliderLabels();
  updateProfileStatus();
}

void NoiseRepellentAudioProcessorEditor::paint(juce::Graphics& g) {
  g.fillAll(juce::Colour(0xff2c313c));

  if (isAdvancedVisible) {
    g.setColour(NoiseRepellentLookAndFeel::kColorPanelBorder);

    float topY = (float)advanced.groupAdvanced.getY() + 22.0f;
    float bottomY = (float)advanced.groupAdvanced.getBottom() - 10.0f;

    // Separator lines between adjacent columns in Advanced Controls (left
    // to right: Reduction Blend, Smoothing, Masking, Whitening,
    // Aggressiveness — low latency swaps Masking for a second pair gap).
    // Bounds are always valid once the panel is laid out, so the midpoints
    // draw unconditionally and the pair list only encodes the adjacency
    // order of the ACTIVE columns (both are visible, no wrap-around).
    std::array<const juce::Slider*, 5> kColumns = {
        &advanced.sliderReductionBlend, &advanced.sliderSmoothing,
        &advanced.sliderMasking, &advanced.sliderWhitening,
        &advanced.sliderAggressiveness};
    int numColumns = lastAdvancedLowLatency ? 4 : 5;
    if (lastAdvancedLowLatency) {
      kColumns[2] = kColumns[3]; // whitening shifts into Masking's slot
      kColumns[3] = kColumns[4]; // aggressiveness into Whitening's slot
    }
    for (int i = 0; i < numColumns - 1; ++i) {
      const float x =
          (float)(kColumns[i]->getRight() + kColumns[i + 1]->getX()) / 2.0f;
      g.drawVerticalLine(juce::roundToInt(x), topY, bottomY);
    }
  }
}

void NoiseRepellentAudioProcessorEditor::paintOverChildren(juce::Graphics& g) {
  if (!wasOfflineRendering)
    return;

  // Darken UI overlay
  auto bounds = getLocalBounds();
  g.setColour(juce::Colour(0xd0161a22));
  g.fillRect(bounds);

  // Center badge card
  constexpr int cardW = 320;
  constexpr int cardH = 96;
  auto cardRect = bounds.withSizeKeepingCentre(cardW, cardH).toFloat();

  g.setColour(juce::Colour(0xee212630));
  g.fillRoundedRectangle(cardRect, 8.0f);

  g.setColour(NoiseRepellentLookAndFeel::kColorPanelBorder);
  g.drawRoundedRectangle(cardRect, 8.0f, 1.5f);

  auto contentArea = cardRect.toNearestInt().reduced(12);

  // Status title
  g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
  g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeBrand,
                              juce::Font::bold));
  g.drawText("OFFLINE RENDERING", contentArea.removeFromTop(24),
             juce::Justification::centred, false);

  contentArea.removeFromTop(4);

  // Subtitle explanation
  g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorFooterText));
  g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                              juce::Font::plain));
  g.drawText("Processing in progress. UI disabled for speed.", contentArea,
             juce::Justification::centred, false);
}

void NoiseRepellentAudioProcessorEditor::resized() {
  // Layout metrics (pixel geometry; values define the look, names the intent)
  constexpr int kOuterMargin = 12;
  constexpr int kHeaderH = 54;
  constexpr int kBrandW = 135;
  constexpr int kBrandH = 26;
  constexpr int kTitleGap = 4;
  constexpr int kPrefsBtn = 22;
  constexpr int kBypassW = 60;
  constexpr int kHeaderBtnH = 24;
  constexpr int kHeaderBtnGap = 6;
  constexpr int kDeltaW = 50;
  constexpr int kAlgoColH = 44;
  constexpr int kAlgoLabelH = 15;
  constexpr int kAlgoLabelGap = 3;
  constexpr int kAlgoComboH = 26;
  constexpr int kGroupPadX = 6;
  constexpr int kGroupPadY = 4;
  constexpr int kGroupTitleH = 10;
  constexpr int kProfileBtnH = 24;
  constexpr int kSectionGap = 8;
  constexpr int kAdvPanelH = 68;
  constexpr int kAdvLabelH = 16;
  constexpr int kAdvLabelGap = 2;
  constexpr int kAdvSliderH = 20;
  constexpr int kInnerPad = 10;
  constexpr int kBankLinkedW = 95;
  constexpr int kBankUnlinkedW = 155;
  constexpr int kBankLabelH = 16;
  constexpr int kBankLabelGap = 2;
  constexpr int kBankBtnH = 22;
  constexpr int kBankBtnGap = 6;
  constexpr int kBankSplitGap = 6;
  constexpr int kCurveResetW = 22;
  constexpr int kCurveResetGap = 2;
  constexpr int kBankGutter = 10;
  constexpr int kAdvTitleH = 12;

  auto area = getLocalBounds().reduced(kOuterMargin);

  // Header (Fixed 54px with spacious headroom for profile controls)
  auto headerArea = area.removeFromTop(kHeaderH);

  // Left Title & Options Block
  header.brandButton.setBounds(
      headerArea.removeFromLeft(kBrandW).withSizeKeepingCentre(kBrandW,
                                                               kBrandH));
  headerArea.removeFromLeft(kTitleGap);
  header.btnPreferences.setBounds(
      headerArea.removeFromLeft(kPrefsBtn).withSizeKeepingCentre(kPrefsBtn,
                                                                 kPrefsBtn));

  // Right Action Buttons Block
  header.btnBypass.setBounds(
      headerArea.removeFromRight(kBypassW).withSizeKeepingCentre(kBypassW,
                                                                 kHeaderBtnH));
  headerArea.removeFromRight(kHeaderBtnGap);
  header.btnDelta.setBounds(
      headerArea.removeFromRight(kDeltaW).withSizeKeepingCentre(kDeltaW,
                                                                kHeaderBtnH));

  constexpr int kHeaderGap = 16;
  headerArea.removeFromLeft(kHeaderGap);

  int availMiddleWidth = headerArea.getWidth() - kHeaderGap;

  constexpr int kLearnW = 72;
  constexpr int kAdaptW = 56;
  constexpr int kArrowW = 18;
  constexpr int kResetW = 22;
  constexpr int kBtnGap = 4;
  int buttonsWidth =
      isAdvancedVisible
          ? (kLearnW + kBtnGap + kAdaptW + kArrowW + kBtnGap + kResetW)
          : (kLearnW + kBtnGap + kAdaptW + kBtnGap + kResetW);

  constexpr int kAlgoWidth = 210;
  int kProfileWidth = buttonsWidth + 12;
  if (int combinedHeaderWidth = kAlgoWidth + kHeaderGap + kProfileWidth;
      availMiddleWidth > combinedHeaderWidth) {
    int leftPad = (availMiddleWidth - combinedHeaderWidth) / 2;
    headerArea.removeFromLeft(leftPad);
  }

  // Processing Engine Column in Header
  auto algoCol = headerArea.removeFromLeft(kAlgoWidth);
  int algoYPad = std::max(0, (algoCol.getHeight() - kAlgoColH) / 2);
  algoCol.removeFromTop(algoYPad);
  header.lblAlgoHeader.setBounds(algoCol.removeFromTop(kAlgoLabelH));
  algoCol.removeFromTop(kAlgoLabelGap);
  header.comboAlgoMode.setBounds(algoCol.removeFromTop(kAlgoComboH));

  headerArea.removeFromLeft(kHeaderGap);

  // Encapsulated Noise Profile Group Box in Header
  auto profileGroupArea = headerArea.removeFromLeft(kProfileWidth);
  profile.groupProfile.setBounds(profileGroupArea);

  auto profileInner = profileGroupArea.reduced(kGroupPadX, kGroupPadY);
  profileInner.removeFromTop(kGroupTitleH);

  auto profileRow = profileInner;

  auto bLearn = profileRow.removeFromLeft(kLearnW);
  profile.btnLearn.setBounds(
      bLearn.withSizeKeepingCentre(kLearnW, kProfileBtnH));
  profileRow.removeFromLeft(kBtnGap);

  auto bAdapt = profileRow.removeFromLeft(kAdaptW);
  profile.btnAdaptiveNoise.setBounds(
      bAdapt.withSizeKeepingCentre(kAdaptW, kProfileBtnH));

  if (isAdvancedVisible) {
    auto bAdaptArrow = profileRow.removeFromLeft(kArrowW);
    profile.btnAdaptiveArrow.setBounds(
        bAdaptArrow.withSizeKeepingCentre(kArrowW, kProfileBtnH));
  }
  profileRow.removeFromLeft(kBtnGap);

  auto bReset = profileRow.removeFromLeft(kResetW);
  profile.btnResetProfile.setBounds(
      bReset.withSizeKeepingCentre(kResetW, kProfileBtnH));

  area.removeFromTop(kSectionGap);

  // Footer bar
  layoutFooterBar(area);

  // Bottom Collapsible Advanced Panel (Reduction Blend, Smoothing, Masking,
  // Whitening, Aggressiveness — low latency shows Reduction Blend instead
  // of Masking)
  if (isAdvancedVisible) {
    auto advArea = area.removeFromBottom(kAdvPanelH);
    advanced.groupAdvanced.setBounds(advArea);

    auto advInner = advArea.reduced(kInnerPad, kGroupPadY);
    advInner.removeFromTop(kAdvTitleH);

    constexpr int kGapW = 14;
    // In low latency the panel swaps the Masking column for the Reduction
    // Blend column (4 columns); otherwise Masking keeps its column (5).
    const int numColumns = lastAdvancedLowLatency ? 4 : 5;
    const int kNumGaps = numColumns - 1;
    // Columns stay full width for the labels, but each slider track is
    // inset so an end-stop thumb never pokes into the separator gap
    // between two adjacent sliders.
    constexpr int kTrackInset = kGapW / 2 + 1;
    int totalAvailWidth = advInner.getWidth() - (kNumGaps * kGapW);
    int itemW = totalAvailWidth / numColumns;

    // 1. Reduction Blend Slider
    auto s1 = advInner.removeFromLeft(itemW);
    advanced.lblReductionBlend.setBounds(s1.removeFromTop(kAdvLabelH));
    s1.removeFromTop(kAdvLabelGap);
    advanced.sliderReductionBlend.setBounds(
        s1.removeFromTop(kAdvSliderH).reduced(kTrackInset, 0));

    // 2. Smoothing Slider
    advInner.removeFromLeft(kGapW);
    auto s2 = advInner.removeFromLeft(itemW);
    advanced.lblSmoothing.setBounds(s2.removeFromTop(kAdvLabelH));
    s2.removeFromTop(kAdvLabelGap);
    advanced.sliderSmoothing.setBounds(
        s2.removeFromTop(kAdvSliderH).reduced(kTrackInset, 0));

    // 3. Masking Slider (hidden in low latency — Reduction Blend takes
    // the panel's first column there)
    if (!lastAdvancedLowLatency) {
      advInner.removeFromLeft(kGapW);
      auto s3 = advInner.removeFromLeft(itemW);
      advanced.lblMasking.setBounds(s3.removeFromTop(kAdvLabelH));
      s3.removeFromTop(kAdvLabelGap);
      advanced.sliderMasking.setBounds(
          s3.removeFromTop(kAdvSliderH).reduced(kTrackInset, 0));
    }

    // 4. Whitening Slider
    advInner.removeFromLeft(kGapW);
    auto s4 = advInner.removeFromLeft(itemW);
    advanced.lblWhitening.setBounds(s4.removeFromTop(kAdvLabelH));
    s4.removeFromTop(kAdvLabelGap);
    advanced.sliderWhitening.setBounds(
        s4.removeFromTop(kAdvSliderH).reduced(kTrackInset, 0));

    // 5. Aggressiveness Slider (fills the remainder, absorbing rounding
    // leftovers from the integer division)
    advInner.removeFromLeft(kGapW);
    auto s5 = advInner;
    advanced.lblAggressiveness.setBounds(s5.removeFromTop(kAdvLabelH));
    s5.removeFromTop(kAdvLabelGap);
    advanced.sliderAggressiveness.setBounds(
        s5.removeFromTop(kAdvSliderH).reduced(kTrackInset, 0));

    area.removeFromBottom(kSectionGap);
  }

  // Main Denoising & Spectrum Canvas
  denoising.groupDenoising.setBounds(area);

  auto denoiseInner = area.reduced(kInnerPad);
  denoiseInner.removeFromTop(kGroupTitleH);

  bool isLinked = !isAdvancedVisible || denoising.btnLink.getToggleState();
  int faderBankWidth = isLinked ? kBankLinkedW : kBankUnlinkedW;

  auto faderBankArea = denoiseInner.removeFromLeft(faderBankWidth);
  denoising.lblReductionHeader.setBounds(
      faderBankArea.removeFromTop(kBankLabelH));
  faderBankArea.removeFromTop(kBankLabelGap);

  if (isAdvancedVisible) {
    denoising.btnLink.setBounds(faderBankArea.removeFromTop(kBankBtnH));
    faderBankArea.removeFromTop(kBankBtnGap);

    // Carve bottom area of faderBankArea for Curve split button [ Curve ][ ↺ ]
    bool isCurveOn = denoising.btnCurveToggle.getToggleState();
    auto faderBottomArea = faderBankArea.removeFromBottom(kBankBtnH);
    faderBankArea.removeFromBottom(kBankBtnGap);

    if (isCurveOn) {
      denoising.btnCurveToggle.setBounds(faderBottomArea.removeFromLeft(
          faderBottomArea.getWidth() - kCurveResetW));
      faderBottomArea.removeFromLeft(kCurveResetGap);
      denoising.btnResetCurve.setBounds(faderBottomArea);
    } else {
      denoising.btnCurveToggle.setBounds(faderBottomArea);
    }
  }

  if (isLinked) {
    denoising.sliderMasterRed.setBounds(faderBankArea);
  } else {
    int colW = (faderBankArea.getWidth() - kBankSplitGap) / 2;

    auto leftCol = faderBankArea.removeFromLeft(colW);
    faderBankArea.removeFromLeft(kBankSplitGap);
    auto rightCol = faderBankArea;

    denoising.lblMasterRed.setBounds(leftCol.removeFromTop(kBankLabelH));
    leftCol.removeFromTop(kBankLabelGap);
    denoising.sliderMasterRed.setBounds(leftCol);

    denoising.lblTonalRed.setBounds(rightCol.removeFromTop(kBankLabelH));
    rightCol.removeFromTop(kBankLabelGap);
    denoising.sliderTonalRed.setBounds(rightCol);
  }

  denoiseInner.removeFromLeft(kBankGutter);

  // Right-side Threshold (Noise Profile Offset) vertical fader bank (symmetric
  // to reduction sliders on left)
  bool isOffsetLinked =
      !isAdvancedVisible || profile.btnLinkOffset.getToggleState();
  int offsetBankWidth = isOffsetLinked ? kBankLinkedW : kBankUnlinkedW;

  auto offsetBankArea = denoiseInner.removeFromRight(offsetBankWidth);
  denoiseInner.removeFromRight(kBankGutter);

  profile.lblOffset.setBounds(offsetBankArea.removeFromTop(kBankLabelH));
  offsetBankArea.removeFromTop(kBankLabelGap);

  if (isAdvancedVisible) {
    profile.btnLinkOffset.setBounds(offsetBankArea.removeFromTop(kBankBtnH));
    offsetBankArea.removeFromTop(kBankBtnGap);
  }

  // Carve bottom area of offsetBankArea for Advanced Controls button below
  // threshold slider
  auto offsetBottomArea = offsetBankArea.removeFromBottom(kBankBtnH);
  offsetBankArea.removeFromBottom(kBankBtnGap);
  header.btnAdvancedToggle.setBounds(offsetBottomArea);

  if (isOffsetLinked) {
    profile.sliderOffset.setBounds(offsetBankArea);
  } else {
    int colW = (offsetBankArea.getWidth() - kBankSplitGap) / 2;

    auto leftCol = offsetBankArea.removeFromLeft(colW);
    offsetBankArea.removeFromLeft(kBankSplitGap);
    auto rightCol = offsetBankArea;

    profile.lblMasterOffset.setBounds(leftCol.removeFromTop(kBankLabelH));
    leftCol.removeFromTop(kBankLabelGap);
    profile.sliderOffset.setBounds(leftCol);

    profile.lblTonalOffset.setBounds(rightCol.removeFromTop(kBankLabelH));
    rightCol.removeFromTop(kBankLabelGap);
    profile.sliderTonalOffset.setBounds(rightCol);
  }

  spectralVisualizer.setBounds(denoiseInner);

  header.btnPreferences.toFront(false);
  footerTooltipLabel.toFront(false);
}

void NoiseRepellentAudioProcessorEditor::layoutFooterBar(
    juce::Rectangle<int>& area) {
  // Layout metrics (pixel geometry; values define the look, names the intent)
  constexpr int kFooterGap = 8;
  constexpr int kFooterSplitPct = 55;
  constexpr int kFooterHeight = 34;

  const auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  const auto* showHudParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_hud");

  bool showTooltips =
      (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f);
  bool showHud = (showHudParam == nullptr || showHudParam->load() > 0.5f);
  bool showFooter = showTooltips || showHud;

  footerTooltipLabel.setVisible(showTooltips);
  profile.lblProfileStatus.setVisible(showHud);

  if (showFooter) {
    auto footerArea = area.removeFromBottom(kFooterHeight);
    area.removeFromBottom(kFooterGap);

    if (showTooltips && showHud) {
      auto leftFooter = footerArea.removeFromLeft(footerArea.getWidth() *
                                                  kFooterSplitPct / 100);
      footerTooltipLabel.setBounds(leftFooter);
      profile.lblProfileStatus.setBounds(footerArea);
      profile.lblProfileStatus.setJustificationType(
          juce::Justification::centredRight);
    } else if (showTooltips) {
      footerTooltipLabel.setBounds(footerArea);
    } else if (showHud) {
      profile.lblProfileStatus.setBounds(footerArea);
      profile.lblProfileStatus.setJustificationType(
          juce::Justification::centredRight);
    }
  }
}
