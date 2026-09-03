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

namespace {

// Shared tooltip/instruction strings (single source; previously copy-pasted)
const juce::String kTipDefaultInstruction =
    "Loop a noise-only segment in your DAW and click Learn Noise\nto "
    "capture a profile, then adjust reduction level.";
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
const juce::String kTipAlgoMode =
    "How the noise reduction smoothing is computed. Standard is fast and\n"
    "light on CPU; Patch-Based analyzes similar patches for higher quality.";
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
const juce::String kTipAdaptiveMethod =
    "Select Adaptive Estimation Method: SPP-MMSE (best for speech & dynamic "
    "noise),\nBrandt (best for steady hiss & fans), or Martin (best for slow "
    "background).";
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
const juce::String kTipCurveToggle = "Enable per-frequency reduction bias curve.";
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
    "control Patch-Based (High Quality) smoothing.";
const juce::String kTipMasking =
    "Adjust psychoacoustic masking threshold\nto protect quiet musical "
    "transients.";
const juce::String kTipWhitening =
    "Spectral whitening factor to equalize\nthe residual noise floor "
    "spectrum.";
const juce::String kTipMorphingUnavailable =
    "Profile Morphing is available in Advanced Controls\nwhen a manual "
    "learned noise profile is captured.";

} // namespace
#include "PluginProcessor.h"

NoiseRepellentAudioProcessorEditor::NoiseRepellentAudioProcessorEditor(
    NoiseRepellentAudioProcessor& p)
    : AudioProcessorEditor(&p), audioProcessor(p), spectralVisualizer(p) {
  setLookAndFeel(&customLookAndFeel);

  // 1. Native Resizable Window Settings
  setResizable(true, true);
  setResizeLimits(840, 480, 1600, 1000);
  setSize(940, 560);

  // Brand Header (accessible button opening the About box with dependency
  // versions)
  brandButton.setName(NoiseRepellentLookAndFeel::kBrandButtonName);
  brandButton.setButtonText("NOISE REPELLENT");
  brandButton.setTooltip(kTipAbout);
  brandButton.setMouseCursor(juce::MouseCursor::PointingHandCursor);
  brandButton.onClick = [this]() { showAboutBox(); };
  addAndMakeVisible(brandButton);

  // Preferences Header Button
  addAndMakeVisible(btnPreferences);
  btnPreferences.setTooltip(kTipPreferences);
  btnPreferences.setColour(juce::TextButton::buttonColourId,
                           juce::Colour(0xff3a4150));
  btnPreferences.onClick = [this]() {
    auto* tooltipParam =
        audioProcessor.getAPVTS().getParameter("show_tooltips");
    auto* hudParam = audioProcessor.getAPVTS().getParameter("show_hud");
    auto* frameSizeParam =
        audioProcessor.getAPVTS().getParameter("frame_size");

    bool tooltipVal =
        (tooltipParam != nullptr && tooltipParam->getValue() > 0.5f);
    bool hudVal = (hudParam != nullptr && hudParam->getValue() > 0.5f);
    int frameSizeIdx = 2;
    if (auto* choice = dynamic_cast<juce::AudioParameterChoice*>(
            frameSizeParam))
      frameSizeIdx = choice->getIndex();

    juce::PopupMenu menu;

    juce::PopupMenu::Item item1("Show Tooltips");
    item1.itemID = 1;
    item1.isTicked = tooltipVal;
    item1.action = [this, tooltipParam, tooltipVal]() {
      if (tooltipParam != nullptr) {
        tooltipParam->beginChangeGesture();
        tooltipParam->setValueNotifyingHost(tooltipVal ? 0.0f : 1.0f);
        tooltipParam->endChangeGesture();
        updateLayout();
      }
    };
    menu.addItem(item1);

    juce::PopupMenu::Item item2("Show HUD Status");
    item2.itemID = 2;
    item2.isTicked = hudVal;
    item2.action = [this, hudParam, hudVal]() {
      if (hudParam != nullptr) {
        hudParam->beginChangeGesture();
        hudParam->setValueNotifyingHost(hudVal ? 0.0f : 1.0f);
        hudParam->endChangeGesture();
        updateLayout();
      }
    };
    menu.addItem(item2);

    menu.addSeparator();
    menu.addSectionHeader("STFT FRAME SIZE (FFT RESOLUTION)");
    static constexpr const char* kFrameSizeNames[4] = {"23 ms", "32 ms",
                                                       "46 ms", "64 ms"};
    for (int i = 0; i < 4; ++i) {
      juce::PopupMenu::Item frameItem(kFrameSizeNames[i]);
      frameItem.itemID = 3 + i;
      frameItem.isTicked = (i == frameSizeIdx);
      frameItem.action = [this, frameSizeParam, i]() {
        if (frameSizeParam != nullptr) {
          frameSizeParam->beginChangeGesture();
          frameSizeParam->setValueNotifyingHost(static_cast<float>(i) / 3.0f);
          frameSizeParam->endChangeGesture();
        }
      };
      menu.addItem(frameItem);
    }

    menu.setLookAndFeel(&getLookAndFeel());
    menu.showMenuAsync(
        juce::PopupMenu::Options().withTargetComponent(&btnPreferences));
  };

  // Algorithm Mode Label & Dropdown
  lblAlgoHeader.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblAlgoHeader.setColour(juce::Label::textColourId,
                          NoiseRepellentLookAndFeel::kColorDenoising);
  lblAlgoHeader.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblAlgoHeader);

  comboAlgoMode.addItemList(
      {"Standard (Fast & Low CPU)", "Patch-Based (High Quality)"}, 1);
  addAndMakeVisible(comboAlgoMode);
  comboAlgoMode.onChange = [this]() { updateLayout(); };

  // Header Action Toggles
  isAdvancedVisible =
      audioProcessor.getAPVTS().getRawParameterValue("show_advanced")->load() >
      0.5f;
  btnAdvancedToggle.setClickingTogglesState(true);
  addAndMakeVisible(btnAdvancedToggle);
  btnAdvancedToggle.setColour(juce::TextButton::buttonColourId,
                              juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  btnAdvancedToggle.setColour(juce::TextButton::buttonOnColourId,
                              NoiseRepellentLookAndFeel::kColorNoiseProfile);
  btnAdvancedToggle.onClick = [this]() {
    isAdvancedVisible = btnAdvancedToggle.getToggleState();
    updateLayout();
  };

  btnDelta.setClickingTogglesState(true);
  btnDelta.setColour(juce::TextButton::buttonColourId,
                     juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  btnDelta.setColour(juce::TextButton::buttonOnColourId,
                     NoiseRepellentLookAndFeel::kColorNoiseProfile);
  addAndMakeVisible(btnDelta);
  btnDelta.onClick = [this]() { updateLayout(); };

  btnBypass.setClickingTogglesState(true);
  btnBypass.setColour(juce::TextButton::buttonColourId,
                      juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  btnBypass.setColour(juce::TextButton::buttonOnColourId,
                      NoiseRepellentLookAndFeel::kColorNoiseProfile);
  addAndMakeVisible(btnBypass);
  btnBypass.onClick = [this]() { updateLayout(); };

  // Beginner Module 1: Noise Profile
  addAndMakeVisible(groupProfile);
  groupProfile.setText("NOISE PROFILE");
  groupProfile.setColour(juce::GroupComponent::outlineColourId,
                         NoiseRepellentLookAndFeel::kColorPanelBorder);
  groupProfile.setColour(juce::GroupComponent::textColourId,
                         NoiseRepellentLookAndFeel::kColorFineTuning);
  groupProfile.setInterceptsMouseClicks(false, false);

  btnLearn.setClickingTogglesState(true);
  addAndMakeVisible(btnLearn);
  btnLearn.setColour(
      juce::TextButton::buttonColourId,
      juce::Colour(NoiseRepellentLookAndFeel::kColorLearnCTA)); // Prominent Studio Crimson Red CTA
  btnLearn.setColour(juce::TextButton::buttonOnColourId,
                     juce::Colour(NoiseRepellentLookAndFeel::kColorLearnActive)); // Active Learning Red
  btnLearn.onClick = [this]() { updateProfileStatus(); };

  btnAdaptiveNoise.setClickingTogglesState(true);
  btnAdaptiveNoise.setButtonText("Adaptive");
  btnAdaptiveNoise.setTooltip(
      kTipAdaptive);
  addAndMakeVisible(btnAdaptiveNoise);
  btnAdaptiveNoise.setColour(juce::TextButton::buttonColourId,
                             juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
  btnAdaptiveNoise.setColour(juce::TextButton::buttonOnColourId,
                             NoiseRepellentLookAndFeel::kColorDenoising);
  btnAdaptiveNoise.onClick = [this]() {
    updateLayout();
    updateProfileStatus();
  };

  addAndMakeVisible(btnAdaptiveArrow);
  btnAdaptiveArrow.setTooltip(kTipAdaptiveMethod);
  btnAdaptiveArrow.setColour(juce::TextButton::buttonColourId,
                             juce::Colour(0xff353b48));
  btnAdaptiveArrow.onClick = [this]() {
    juce::PopupMenu menu;
    int currentMethod = static_cast<int>(comboMethod.getSelectedId());

    menu.addSectionHeader("ADAPTIVE ESTIMATION METHOD");
    menu.addItem(1, "SPP-MMSE (Unbiased)", true, currentMethod == 1);
    menu.addItem(2, "Brandt (Trimmed Mean)", true, currentMethod == 2);
    menu.addItem(3, "Martin (Minimum Statistics)", true, currentMethod == 3);

    menu.setLookAndFeel(&getLookAndFeel());
    menu.showMenuAsync(
        juce::PopupMenu::Options().withTargetComponent(&btnAdaptiveArrow),
        [this](int result) {
          if (result >= 1 && result <= 3) {
            int methodIdx = result - 1;
            comboMethod.setSelectedId(result, juce::sendNotification);
            if (auto* p =
                    audioProcessor.getAPVTS().getParameter("adaptive_method"))
              p->setValueNotifyingHost(static_cast<float>(methodIdx) / 2.0f);
            if (auto* p =
                    audioProcessor.getAPVTS().getParameter("adaptive_noise"))
              p->setValueNotifyingHost(1.0f);
            updateLayout();
            updateProfileStatus();
          }
        });
  };

  addAndMakeVisible(btnResetProfile);
  btnResetProfile.setTooltip(kTipResetProfile);
  btnResetProfile.onClick = [this]() {
    audioProcessor.resetNoiseProfile();
    if (auto* pLearn = audioProcessor.getAPVTS().getParameter("learn_noise"))
      pLearn->setValueNotifyingHost(0.0f);
    if (auto* pAdapt = audioProcessor.getAPVTS().getParameter("adaptive_noise"))
      pAdapt->setValueNotifyingHost(0.0f);
    updateLayout();
    updateProfileStatus();
  };

  // Profile Status Label (HUD Overlay on Chart)
  addAndMakeVisible(lblProfileStatus);
  lblProfileStatus.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblProfileStatus.setColour(juce::Label::textColourId,
                             NoiseRepellentLookAndFeel::kColorNoiseProfile);
  lblProfileStatus.setJustificationType(juce::Justification::centredRight);
  lblProfileStatus.setText("STATUS: STATIONARY NOISE PROFILE",
                           juce::dontSendNotification);

  // Beginner Module 2: Denoising & Resizable Canvas
  addAndMakeVisible(groupDenoising);
  groupDenoising.setText("");
  groupDenoising.setColour(juce::GroupComponent::outlineColourId,
                           NoiseRepellentLookAndFeel::kColorPanelBorder);
  groupDenoising.setColour(juce::GroupComponent::textColourId,
                           NoiseRepellentLookAndFeel::kColorFineTuning);
  groupDenoising.setInterceptsMouseClicks(false, false);

  lblReductionHeader.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblReductionHeader.setColour(juce::Label::textColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
  lblReductionHeader.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblReductionHeader);

  addAndMakeVisible(btnLink);
  // Attachment must be created before onClick is assigned: ButtonAttachment
  // calls setToggleState(sendNotification) during construction to sync from
  // APVTS, which would fire onClick and reset reduction params to default.
  attachLink = std::make_unique<ButtonAttachment>(audioProcessor.getAPVTS(),
                                                  "link_reduction", btnLink);
  btnLink.setButtonText(btnLink.getToggleState() ? "Linked" : "Unlinked");
  btnLink.onClick = [this]() {
    bool isLinked = btnLink.getToggleState();
    btnLink.setButtonText(isLinked ? "Linked" : "Unlinked");
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
  };

  btnCurveToggle.setClickingTogglesState(true);
  btnCurveToggle.setTooltip(kTipCurveToggle);
  addAndMakeVisible(btnCurveToggle);
  btnCurveToggle.onClick = [this]() { updateLayout(); };

  btnResetCurve.setTooltip(kTipResetCurve);
  addAndMakeVisible(btnResetCurve);
  btnResetCurve.onClick = [this]() {
    audioProcessor.resetCurveNodes();
    spectralVisualizer.repaint();
  };

  lblMasterRed.setText("BROADBAND", juce::dontSendNotification);
  lblMasterRed.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblMasterRed.setColour(juce::Label::textColourId,
                         NoiseRepellentLookAndFeel::kColorDenoising);
  lblMasterRed.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblMasterRed);

  lblTonalRed.setText("TONAL", juce::dontSendNotification);
  lblTonalRed.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblTonalRed.setColour(juce::Label::textColourId,
                        NoiseRepellentLookAndFeel::kColorDenoising);
  lblTonalRed.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblTonalRed);

  sliderMasterRed.setColour(juce::Slider::rotarySliderFillColourId,
                            NoiseRepellentLookAndFeel::kColorDenoising);
  sliderTonalRed.setColour(juce::Slider::rotarySliderFillColourId,
                           NoiseRepellentLookAndFeel::kColorDenoising);

  addAndMakeVisible(sliderMasterRed);
  addAndMakeVisible(sliderTonalRed);

  addAndMakeVisible(spectralVisualizer);

  // Collapsible Module 3: Advanced Controls
  addAndMakeVisible(groupAdvanced);
  groupAdvanced.setText("ADVANCED CONTROLS");
  groupAdvanced.setColour(juce::GroupComponent::outlineColourId,
                          NoiseRepellentLookAndFeel::kColorPanelBorder);
  groupAdvanced.setColour(juce::GroupComponent::textColourId,
                          NoiseRepellentLookAndFeel::kColorFineTuning);
  groupAdvanced.setInterceptsMouseClicks(false, false);

  comboMethod.addItemList({"SPP-MMSE (Unbiased)", "Brandt (Trimmed Mean)",
                           "Martin (Min Statistics)"},
                          1);
  addAndMakeVisible(comboMethod);

  lblMethod.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                                      juce::Font::bold));
  lblMethod.setColour(juce::Label::textColourId,
                      NoiseRepellentLookAndFeel::kColorDenoising);
  lblMethod.setJustificationType(juce::Justification::left);
  addAndMakeVisible(lblMethod);

  lblOffset.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                                      juce::Font::bold));
  lblOffset.setColour(juce::Label::textColourId,
                      NoiseRepellentLookAndFeel::kColorDenoising);
  lblOffset.setJustificationType(juce::Justification::centred);
  sliderOffset.setColour(juce::Slider::rotarySliderFillColourId,
                         NoiseRepellentLookAndFeel::kColorNoiseProfile);
  addAndMakeVisible(sliderOffset);
  addAndMakeVisible(lblOffset);

  lblMasterOffset.setText("BROADBAND", juce::dontSendNotification);
  lblMasterOffset.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblMasterOffset.setColour(juce::Label::textColourId,
                            NoiseRepellentLookAndFeel::kColorNoiseProfile);
  lblMasterOffset.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblMasterOffset);

  lblTonalOffset.setText("TONAL", juce::dontSendNotification);
  lblTonalOffset.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblTonalOffset.setColour(juce::Label::textColourId,
                           NoiseRepellentLookAndFeel::kColorTonalPeaks);
  lblTonalOffset.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblTonalOffset);

  sliderTonalOffset.setColour(juce::Slider::rotarySliderFillColourId,
                              NoiseRepellentLookAndFeel::kColorTonalPeaks);

  addAndMakeVisible(sliderTonalOffset);

  addAndMakeVisible(btnLinkOffset);
  // Attachment must be created before onClick is assigned (same reasoning as
  // the reduction link button above).
  attachLinkOffset = std::make_unique<ButtonAttachment>(
      audioProcessor.getAPVTS(), "link_threshold_offset", btnLinkOffset);
  btnLinkOffset.setButtonText(btnLinkOffset.getToggleState() ? "Linked"
                                                             : "Unlinked");
  btnLinkOffset.onClick = [this]() {
    bool isLinked = btnLinkOffset.getToggleState();
    btnLinkOffset.setButtonText(isLinked ? "Linked" : "Unlinked");
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
  };

  lblAggressiveness.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblAggressiveness.setColour(juce::Label::textColourId,
                              NoiseRepellentLookAndFeel::kColorDenoising);
  lblAggressiveness.setJustificationType(juce::Justification::centred);
  sliderAggressiveness.setColour(juce::Slider::rotarySliderFillColourId,
                                 NoiseRepellentLookAndFeel::kColorDenoising);
  addAndMakeVisible(sliderAggressiveness);
  addAndMakeVisible(lblAggressiveness);

  sliderSmoothing.setColour(juce::Slider::rotarySliderFillColourId,
                            NoiseRepellentLookAndFeel::kColorDenoising);
  sliderMasking.setColour(juce::Slider::rotarySliderFillColourId,
                          NoiseRepellentLookAndFeel::kColorDenoising);
  sliderWhitening.setColour(juce::Slider::rotarySliderFillColourId,
                            NoiseRepellentLookAndFeel::kColorDenoising);

  addAndMakeVisible(sliderSmoothing);
  addAndMakeVisible(sliderMasking);
  addAndMakeVisible(sliderWhitening);

  lblSmoothing.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblSmoothing.setColour(juce::Label::textColourId,
                         NoiseRepellentLookAndFeel::kColorDenoising);
  lblSmoothing.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblSmoothing);

  lblMasking.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblMasking.setColour(juce::Label::textColourId,
                       NoiseRepellentLookAndFeel::kColorDenoising);
  lblMasking.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblMasking);

  lblWhitening.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblWhitening.setColour(juce::Label::textColourId,
                         NoiseRepellentLookAndFeel::kColorDenoising);
  lblWhitening.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblWhitening);

  // Advanced Controls Toggle Button
  btnAdvancedToggle.setClickingTogglesState(true);
  btnAdvancedToggle.setButtonText("ADVANCED");
  btnAdvancedToggle.setTooltip(kTipAdvancedToggle);
  addAndMakeVisible(btnAdvancedToggle);
  btnAdvancedToggle.onClick = [this]() {
    isAdvancedVisible = btnAdvancedToggle.getToggleState();
    updateLayout();
  };

  // Parameter APVTS Attachments
  auto& apvts = audioProcessor.getAPVTS();

  attachAlgoMode = std::make_unique<ComboBoxAttachment>(apvts, "algorithm_mode",
                                                        comboAlgoMode);
  attachLearn =
      std::make_unique<ButtonAttachment>(apvts, "learn_noise", btnLearn);
  attachAdaptive = std::make_unique<ButtonAttachment>(apvts, "adaptive_noise",
                                                      btnAdaptiveNoise);
  attachCurveToggle = std::make_unique<ButtonAttachment>(
      apvts, "reduction_curve_enabled", btnCurveToggle);
  attachDelta =
      std::make_unique<ButtonAttachment>(apvts, "residual_listen", btnDelta);
  attachBypass = std::make_unique<ButtonAttachment>(apvts, "bypass", btnBypass);
  attachShowAdvanced = std::make_unique<ButtonAttachment>(
      apvts, "show_advanced", btnAdvancedToggle);

  attachMethod = std::make_unique<ComboBoxAttachment>(apvts, "adaptive_method",
                                                      comboMethod);

  attachMasterRed = std::make_unique<SliderAttachment>(
      apvts, "reduction_amount", sliderMasterRed);
  attachTonalRed = std::make_unique<SliderAttachment>(apvts, "tonal_reduction",
                                                      sliderTonalRed);
  attachOffset = std::make_unique<SliderAttachment>(
      apvts, "noise_profile_offset", sliderOffset);
  attachTonalOffset = std::make_unique<SliderAttachment>(
      apvts, "tonal_noise_profile_offset", sliderTonalOffset);
  attachSmoothing = std::make_unique<SliderAttachment>(
      apvts, "smoothing_factor", sliderSmoothing);
  attachMasking =
      std::make_unique<SliderAttachment>(apvts, "masking_depth", sliderMasking);
  attachWhitening = std::make_unique<SliderAttachment>(
      apvts, "whitening_factor", sliderWhitening);
  attachAggressiveness = std::make_unique<SliderAttachment>(
      apvts, "aggressiveness", sliderAggressiveness);

  // Control Tooltip Descriptions
  btnPreferences.setTooltip(kTipPreferencesMenu);
  comboAlgoMode.setTooltip(
kTipAlgoMode);
  lblAlgoHeader.setTooltip(
kTipAlgoMode);
  btnAdvancedToggle.setTooltip(kTipAdvancedShow);
  btnLearn.setTooltip(kTipLearnInitial);
  btnResetProfile.setTooltip(kTipResetProfileClear);
  btnAdaptiveNoise.setTooltip(
      kTipAdaptive);
  btnAdaptiveArrow.setTooltip(kTipAdaptiveMethod);
  comboMethod.setTooltip(kTipAdaptiveMethod);
  lblMethod.setTooltip(kTipAdaptiveMethod);

  sliderMasterRed.setTooltip(kTipMasterReduction);
  lblMasterRed.setTooltip(kTipMasterReduction);
  lblReductionHeader.setTooltip(kTipMasterReduction);
  sliderTonalRed.setTooltip(kTipTonalReduction);
  lblTonalRed.setTooltip(kTipTonalReduction);

  sliderOffset.setTooltip(kTipThreshold);
  lblOffset.setTooltip(kTipThreshold);
  sliderAggressiveness.setTooltip(kTipAggressiveness);
  lblAggressiveness.setTooltip(kTipAggressiveness);
  btnDelta.setTooltip(kTipDelta);
  btnBypass.setTooltip(kTipBypass);

  sliderSmoothing.setTooltip(kTipSmoothing);
  sliderMasking.setTooltip(kTipMasking);
  sliderWhitening.setTooltip(kTipWhitening);
  btnLink.setTooltip(kTipLink);
  btnCurveToggle.setTooltip(kTipCurveOverlay);
  btnResetCurve.setTooltip(kTipResetCurve);

  // Initial Layout update
  isAdvancedVisible = btnAdvancedToggle.getToggleState();
  updateLayout();
  updateSliderLabels();
  startTimerHz(30);

  // Initial default instruction text
  footerTooltipLabel.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeTooltip, juce::Font::plain));
  footerTooltipLabel.setColour(juce::Label::textColourId,
                               juce::Colour(NoiseRepellentLookAndFeel::kColorFooterText));
  footerTooltipLabel.setJustificationType(juce::Justification::left);
  addAndMakeVisible(footerTooltipLabel);

  auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  if (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f) {
    footerTooltipLabel.setText(
        kTipDefaultInstruction,
        juce::dontSendNotification);
  }

  // Register parameter listener for show_tooltips
  audioProcessor.getAPVTS().addParameterListener("show_tooltips", this);

  // Register mouse listener AFTER all child components are added
  addMouseListener(this, true);
}

void NoiseRepellentAudioProcessorEditor::mouseEnter(
    const juce::MouseEvent& event) {
  auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  if (showTooltipsParam != nullptr && showTooltipsParam->load() < 0.5f) {
    footerTooltipLabel.setText({}, juce::dontSendNotification);
    return;
  }

  juce::Component* comp = event.originalComponent;
  if (comp == nullptr || comp == this || comp == &groupProfile ||
      comp == &groupDenoising || comp == &groupAdvanced) {
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
  footerTooltipLabel.setText(
      kTipDefaultInstruction,
      juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::showAboutBox() {
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
          << (runtime != nullptr ? runtime : "unknown");
  if (runtime != nullptr && juce::String(runtime) != expectedBanner)
    message << "\n\nNote: the loaded library differs from the build version.";
  juce::AlertWindow::showMessageBoxAsync(juce::MessageBoxIconType::InfoIcon,
                                         "About Noise Repellent", message, "OK");
}

void NoiseRepellentAudioProcessorEditor::mouseMove(
    const juce::MouseEvent& event) {
  mouseEnter(event);
}

void NoiseRepellentAudioProcessorEditor::mouseExit(
    const juce::MouseEvent& /*event*/) {
  auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  if (showTooltipsParam != nullptr && showTooltipsParam->load() > 0.5f) {
    footerTooltipLabel.setText(
        kTipDefaultInstruction,
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
    tooltipStateDirty.store(true, std::memory_order_relaxed);
  }
}

void NoiseRepellentAudioProcessorEditor::handleAsyncUpdate() {
  auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  if (showTooltipsParam != nullptr && showTooltipsParam->load() > 0.5f) {
    footerTooltipLabel.setText(
        kTipDefaultInstruction,
        juce::dontSendNotification);
  } else {
    footerTooltipLabel.setText({}, juce::dontSendNotification);
  }
  updateLayout();
}

void NoiseRepellentAudioProcessorEditor::updateLayout() {
  bool isLinked = btnLink.getToggleState();
  bool isCurveOn = btnCurveToggle.getToggleState();

  sliderMasterRed.setVisible(true);
  sliderTonalRed.setVisible(isAdvancedVisible && !isLinked);

  lblMasterRed.setVisible(isAdvancedVisible && !isLinked);
  lblTonalRed.setVisible(isAdvancedVisible && !isLinked);

  btnLink.setVisible(isAdvancedVisible);
  btnCurveToggle.setVisible(isAdvancedVisible);
  btnResetCurve.setVisible(isAdvancedVisible && isCurveOn);

  sliderOffset.setVisible(true);
  lblOffset.setVisible(true);

  bool isOffsetLinked = !isAdvancedVisible || btnLinkOffset.getToggleState();
  sliderTonalOffset.setVisible(isAdvancedVisible && !isOffsetLinked);
  lblMasterOffset.setVisible(isAdvancedVisible && !isOffsetLinked);
  lblTonalOffset.setVisible(isAdvancedVisible && !isOffsetLinked);
  btnLinkOffset.setVisible(isAdvancedVisible);

  sliderAggressiveness.setVisible(isAdvancedVisible);
  lblAggressiveness.setVisible(isAdvancedVisible);

  groupAdvanced.setVisible(isAdvancedVisible);
  btnLearn.setVisible(true);
  btnAdaptiveNoise.setVisible(true);
  btnAdaptiveArrow.setVisible(isAdvancedVisible);
  btnResetProfile.setVisible(true);
  comboMethod.setVisible(false);
  lblMethod.setVisible(false);
  sliderSmoothing.setVisible(isAdvancedVisible);
  lblSmoothing.setVisible(isAdvancedVisible);
  sliderMasking.setVisible(isAdvancedVisible);
  lblMasking.setVisible(isAdvancedVisible);
  sliderWhitening.setVisible(isAdvancedVisible);
  lblWhitening.setVisible(isAdvancedVisible);

  spectralVisualizer.setAdvancedControlsVisible(isAdvancedVisible);

  btnAdvancedToggle.setButtonText(juce::CharPointer_UTF8(
      isAdvancedVisible ? "ADVANCED \xe2\x96\xb2" : "ADVANCED \xe2\x96\xbc"));

  auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  auto* showHudParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_hud");

  bool showTooltips =
      (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f);
  bool showHud = (showHudParam == nullptr || showHudParam->load() > 0.5f);

  footerTooltipLabel.setVisible(showTooltips);
  lblProfileStatus.setVisible(showHud);

  resized();
}

void NoiseRepellentAudioProcessorEditor::updateSliderLabels() {
  bool isCompact = getWidth() < 900;
  lblOffset.setText(isCompact ? "THRESH" : "THRESHOLD",
                    juce::dontSendNotification);
  lblAggressiveness.setText(
      "AGGRESSIVENESS: " + juce::String(sliderAggressiveness.getValue(), 1),
      juce::dontSendNotification);
  lblSmoothing.setText(
      "SMOOTHING: " +
          juce::String(static_cast<int>(sliderSmoothing.getValue())) + "%",
      juce::dontSendNotification);
  lblMasking.setText(
      isCompact
          ? "MASKING: " +
                juce::String(static_cast<int>(sliderMasking.getValue())) + "%"
          : "MASKING PROTECT: " +
                juce::String(static_cast<int>(sliderMasking.getValue())) + "%",
      juce::dontSendNotification);
  lblWhitening.setText(
      "WHITENING: " +
          juce::String(static_cast<int>(sliderWhitening.getValue())) + "%",
      juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::updateProfileStatus() {
  auto* adaptParam =
      audioProcessor.getAPVTS().getRawParameterValue("adaptive_noise");
  bool isAdaptive = adaptParam != nullptr && adaptParam->load() > 0.5f;
  bool isLearning = btnLearn.getToggleState();
  bool hasProfile = audioProcessor.hasNoiseProfile();

  if (isLearning) {
    btnLearn.setButtonText("Learning...");
    btnLearn.setColour(juce::TextButton::buttonColourId,
                       juce::Colour(NoiseRepellentLookAndFeel::kColorLearnActive)); // Active Learning Red
    btnLearn.setColour(juce::TextButton::buttonOnColourId,
                       juce::Colour(NoiseRepellentLookAndFeel::kColorLearnActive));
    btnLearn.setTooltip(kTipLearnCapturing);
  } else if (hasProfile) {
    btnLearn.setButtonText("+ Learn");
    btnLearn.setColour(juce::TextButton::buttonColourId,
                       NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btnLearn.setColour(juce::TextButton::buttonOnColourId,
                       NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btnLearn.setTooltip(kTipLearnAccumulate);
  } else {
    btnLearn.setButtonText("Learn Noise");
    btnLearn.setColour(juce::TextButton::buttonColourId,
                       juce::Colour(NoiseRepellentLookAndFeel::kColorLearnCTA)); // Prominent Crimson Red CTA
    btnLearn.setColour(juce::TextButton::buttonOnColourId,
                       juce::Colour(NoiseRepellentLookAndFeel::kColorLearnActive)); // Active Learning Red
    btnLearn.setTooltip(kTipLearnDefault);
  }
  btnLearn.repaint();

  if (isAdaptive) {
    btnAdaptiveNoise.setColour(juce::TextButton::buttonColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
    btnAdaptiveNoise.setColour(juce::TextButton::buttonOnColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
  } else {
    btnAdaptiveNoise.setColour(juce::TextButton::buttonColourId,
                               juce::Colour(NoiseRepellentLookAndFeel::kColorButtonOff));
    btnAdaptiveNoise.setColour(juce::TextButton::buttonOnColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
  }
  btnAdaptiveNoise.repaint();

  // When bypassed, disable everything except the Bypass button itself
  bool isBypassed = btnBypass.getToggleState();
  bool pluginActive = !isBypassed;

  // Header controls
  btnDelta.setEnabled(pluginActive);
  btnAdvancedToggle.setEnabled(pluginActive);
  comboAlgoMode.setEnabled(pluginActive);
  lblAlgoHeader.setEnabled(pluginActive);

  // Noise Profile box: enabled whenever plugin is active
  bool profileEnabled = pluginActive;
  groupProfile.setEnabled(pluginActive);
  btnLearn.setEnabled(profileEnabled);
  btnAdaptiveNoise.setEnabled(profileEnabled);
  btnAdaptiveArrow.setEnabled(profileEnabled);
  btnResetProfile.setEnabled(profileEnabled && (hasProfile || isAdaptive));
  sliderOffset.setEnabled(pluginActive);
  lblOffset.setEnabled(pluginActive);

  // Aggressiveness (Profile Morphing) enabled in Noise Profile box when a
  // manual profile exists and Advanced Controls is ON
  bool aggressivenessEnabled = pluginActive && isAdvancedVisible && hasProfile;
  sliderAggressiveness.setEnabled(aggressivenessEnabled);
  lblAggressiveness.setEnabled(aggressivenessEnabled);

  if (aggressivenessEnabled) {
    sliderAggressiveness.setTooltip(kTipAggressiveness);
    lblAggressiveness.setTooltip(kTipAggressiveness);
  } else {
    sliderAggressiveness.setTooltip(kTipMorphingUnavailable);
    lblAggressiveness.setTooltip(kTipMorphingUnavailable);
  }

  // Reduction controls
  bool canUnlink = (!isAdaptive || hasProfile);
  bool allowUnlink = pluginActive && canUnlink;
  btnLink.setEnabled(allowUnlink);

  if (!canUnlink && !btnLink.getToggleState()) {
    btnLink.setToggleState(true, juce::dontSendNotification);
    if (auto* linkParam =
            audioProcessor.getAPVTS().getParameter("link_reduction"))
      linkParam->setValueNotifyingHost(1.0f);
  }

  if (canUnlink) {
    btnLink.setTooltip(kTipLink);
  } else {
    btnLink.setTooltip(kTipLinkUnlinkedDisabled);
  }

  lblReductionHeader.setEnabled(pluginActive);
  sliderMasterRed.setEnabled(pluginActive);
  lblMasterRed.setEnabled(pluginActive);

  bool isLinked = btnLink.getToggleState();
  bool tonalEnabled = pluginActive && !isLinked && allowUnlink;
  sliderTonalRed.setEnabled(tonalEnabled);
  lblTonalRed.setEnabled(tonalEnabled);

  const juce::String kTipMasterReduction =
      isLinked
          ? "Adjust noise reduction level in decibels\n(0 to 40 dB across all "
            "bands)."
          : "Adjust broadband noise reduction level in decibels\n(0 to 40 dB).";
  const juce::String kTipTonalReduction =
      "Adjust reduction level for tonal noise components\n(0 to 40 dB).";

  lblReductionHeader.setTooltip(kTipMasterReduction);
  sliderMasterRed.setTooltip(kTipMasterReduction);
  lblMasterRed.setTooltip(kTipMasterReduction);
  sliderTonalRed.setTooltip(kTipTonalReduction);
  lblTonalRed.setTooltip(kTipTonalReduction);

  // Threshold Offset controls
  btnLinkOffset.setEnabled(allowUnlink);

  if (!canUnlink && !btnLinkOffset.getToggleState()) {
    btnLinkOffset.setToggleState(true, juce::dontSendNotification);
    btnLinkOffset.setButtonText("Linked");
    if (auto* linkOffsetParam =
            audioProcessor.getAPVTS().getParameter("link_threshold_offset"))
      linkOffsetParam->setValueNotifyingHost(1.0f);
  }

  bool isOffsetLinked = btnLinkOffset.getToggleState();
  bool tonalOffsetEnabled = pluginActive && !isOffsetLinked && allowUnlink;
  sliderTonalOffset.setEnabled(tonalOffsetEnabled);
  lblTonalOffset.setEnabled(tonalOffsetEnabled);


  lblOffset.setTooltip((isOffsetLinked ? kTipThreshold : kTipMasterOffsetUnlinked));
  sliderOffset.setTooltip((isOffsetLinked ? kTipThreshold : kTipMasterOffsetUnlinked));
  lblMasterOffset.setTooltip((isOffsetLinked ? kTipThreshold : kTipMasterOffsetUnlinked));
  sliderTonalOffset.setTooltip(kTipTonalOffset);
  lblTonalOffset.setTooltip(kTipTonalOffset);

  btnCurveToggle.setEnabled(pluginActive);
  btnResetCurve.setEnabled(pluginActive && btnCurveToggle.getToggleState());

  // Advanced controls
  sliderSmoothing.setEnabled(pluginActive);
  sliderMasking.setEnabled(pluginActive);
  sliderWhitening.setEnabled(pluginActive);
  lblSmoothing.setEnabled(pluginActive);
  lblMasking.setEnabled(pluginActive);
  lblWhitening.setEnabled(pluginActive);
  comboMethod.setEnabled(pluginActive);
  lblMethod.setEnabled(pluginActive);
  groupAdvanced.setEnabled(pluginActive);

  bool is2D = (comboAlgoMode.getSelectedItemIndex() == 1);
  const juce::String smoothingTip = is2D ? kTipSmoothing2D : kTipSmoothing;

  sliderSmoothing.setTooltip(smoothingTip);
  lblSmoothing.setTooltip(smoothingTip);

  // Status label (HUD overlay on FFT spectrum chart)
  bool isDelta = btnDelta.getToggleState();

  if (isBypassed) {
    lblProfileStatus.setText("STATUS: BYPASSED", juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               juce::Colour(NoiseRepellentLookAndFeel::kColorInactiveText));
  } else if (isDelta) {
    lblProfileStatus.setText("STATUS: DELTA (RESIDUAL LISTEN)",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               NoiseRepellentLookAndFeel::kColorTonalPeaks);
  } else if (isLearning) {
    lblProfileStatus.setText("STATUS: CAPTURING NOISE PROFILE...",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               juce::Colour(NoiseRepellentLookAndFeel::kColorLearnActive));
  } else if (isAdaptive && hasProfile) {
    lblProfileStatus.setText("STATUS: HYBRID (ADAPTIVE + PROFILE)",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               NoiseRepellentLookAndFeel::kColorNoiseProfile);
  } else if (isAdaptive) {
    lblProfileStatus.setText("STATUS: CONTINUOUS ADAPTIVE",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
  } else if (hasProfile) {
    lblProfileStatus.setText("STATUS: STATIONARY NOISE PROFILE",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               NoiseRepellentLookAndFeel::kColorNoiseProfile);
  } else {
    lblProfileStatus.setText("STATUS: NO PROFILE (PASS-THROUGH)",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               juce::Colour(NoiseRepellentLookAndFeel::kColorInactiveText));
  }
}

void NoiseRepellentAudioProcessorEditor::timerCallback() {
  if (tooltipStateDirty.exchange(false, std::memory_order_relaxed)) {
    triggerAsyncUpdate();
  }

  bool isOffline = audioProcessor.isNonRealtime();
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

    float topY = (float)groupAdvanced.getY() + 22.0f;
    float bottomY = (float)groupAdvanced.getBottom() - 10.0f;

    // 3 Separator lines between 4 controls in Advanced Controls
    if (sliderSmoothing.getWidth() > 0 && sliderMasking.getWidth() > 0) {
      float x0 =
          (float)(sliderSmoothing.getRight() + sliderMasking.getX()) / 2.0f;
      g.drawVerticalLine(juce::roundToInt(x0), topY, bottomY);
    }
    if (sliderMasking.getWidth() > 0 && sliderWhitening.getWidth() > 0) {
      float x1 =
          (float)(sliderMasking.getRight() + sliderWhitening.getX()) / 2.0f;
      g.drawVerticalLine(juce::roundToInt(x1), topY, bottomY);
    }
    if (sliderWhitening.getWidth() > 0 && sliderAggressiveness.getWidth() > 0) {
      float x2 =
          (float)(sliderWhitening.getRight() + sliderAggressiveness.getX()) /
          2.0f;
      g.drawVerticalLine(juce::roundToInt(x2), topY, bottomY);
    }
  }
}

void NoiseRepellentAudioProcessorEditor::paintOverChildren(juce::Graphics& g) {
  if (!audioProcessor.isNonRealtime())
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
  constexpr int kFooterGap = 8;
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
  constexpr int kFooterSplitPct = 55;
  constexpr int kAdvTitleH = 12;

  auto area = getLocalBounds().reduced(kOuterMargin);

  // Header (Fixed 54px with spacious headroom for profile controls)
  auto headerArea = area.removeFromTop(kHeaderH);

  // Left Title & Options Block
  brandButton.setBounds(
      headerArea.removeFromLeft(kBrandW).withSizeKeepingCentre(kBrandW, kBrandH));
  headerArea.removeFromLeft(kTitleGap);
  btnPreferences.setBounds(
      headerArea.removeFromLeft(kPrefsBtn).withSizeKeepingCentre(kPrefsBtn, kPrefsBtn));

  // Right Action Buttons Block
  btnBypass.setBounds(
      headerArea.removeFromRight(kBypassW).withSizeKeepingCentre(kBypassW, kHeaderBtnH));
  headerArea.removeFromRight(kHeaderBtnGap);
  btnDelta.setBounds(
      headerArea.removeFromRight(kDeltaW).withSizeKeepingCentre(kDeltaW, kHeaderBtnH));

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
  int combinedHeaderWidth = kAlgoWidth + kHeaderGap + kProfileWidth;

  if (availMiddleWidth > combinedHeaderWidth) {
    int leftPad = (availMiddleWidth - combinedHeaderWidth) / 2;
    headerArea.removeFromLeft(leftPad);
  }

  // Processing Engine Column in Header
  auto algoCol = headerArea.removeFromLeft(kAlgoWidth);
  int algoYPad = std::max(0, (algoCol.getHeight() - kAlgoColH) / 2);
  algoCol.removeFromTop(algoYPad);
  lblAlgoHeader.setBounds(algoCol.removeFromTop(kAlgoLabelH));
  algoCol.removeFromTop(kAlgoLabelGap);
  comboAlgoMode.setBounds(algoCol.removeFromTop(kAlgoComboH));

  headerArea.removeFromLeft(kHeaderGap);

  // Encapsulated Noise Profile Group Box in Header
  auto profileGroupArea = headerArea.removeFromLeft(kProfileWidth);
  groupProfile.setBounds(profileGroupArea);

  auto profileInner = profileGroupArea.reduced(kGroupPadX, kGroupPadY);
  profileInner.removeFromTop(kGroupTitleH);

  auto profileRow = profileInner;

  auto bLearn = profileRow.removeFromLeft(kLearnW);
  btnLearn.setBounds(bLearn.withSizeKeepingCentre(kLearnW, kProfileBtnH));
  profileRow.removeFromLeft(kBtnGap);

  auto bAdapt = profileRow.removeFromLeft(kAdaptW);
  btnAdaptiveNoise.setBounds(bAdapt.withSizeKeepingCentre(kAdaptW, kProfileBtnH));

  if (isAdvancedVisible) {
    auto bAdaptArrow = profileRow.removeFromLeft(kArrowW);
    btnAdaptiveArrow.setBounds(bAdaptArrow.withSizeKeepingCentre(kArrowW, kProfileBtnH));
  }
  profileRow.removeFromLeft(kBtnGap);

  auto bReset = profileRow.removeFromLeft(kResetW);
  btnResetProfile.setBounds(bReset.withSizeKeepingCentre(kResetW, kProfileBtnH));

  area.removeFromTop(kSectionGap);

  // Footer bar
  auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  auto* showHudParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_hud");

  bool showTooltips =
      (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f);
  bool showHud = (showHudParam == nullptr || showHudParam->load() > 0.5f);
  bool showFooter = showTooltips || showHud;

  footerTooltipLabel.setVisible(showTooltips);
  lblProfileStatus.setVisible(showHud);

  if (showFooter) {
    constexpr int kFooterHeight = 34;
    auto footerArea = area.removeFromBottom(kFooterHeight);
    area.removeFromBottom(kFooterGap);

    if (showTooltips && showHud) {
      auto leftFooter =
          footerArea.removeFromLeft(footerArea.getWidth() * kFooterSplitPct / 100);
      footerTooltipLabel.setBounds(leftFooter);
      lblProfileStatus.setBounds(footerArea);
      lblProfileStatus.setJustificationType(juce::Justification::centredRight);
    } else if (showTooltips) {
      footerTooltipLabel.setBounds(footerArea);
    } else if (showHud) {
      lblProfileStatus.setBounds(footerArea);
      lblProfileStatus.setJustificationType(juce::Justification::centredRight);
    }
  }

  // Bottom Collapsible Advanced Panel (4 Controls: Smoothing, Masking,
  // Whitening, Aggressiveness)
  if (isAdvancedVisible) {
    auto advArea = area.removeFromBottom(kAdvPanelH);
    groupAdvanced.setBounds(advArea);

    auto advInner = advArea.reduced(kInnerPad, kGroupPadY);
    advInner.removeFromTop(kAdvTitleH);

    constexpr int kNumGaps = 3;
    constexpr int kGapW = 14;
    int totalAvailWidth = advInner.getWidth() - (kNumGaps * kGapW);
    int itemW = totalAvailWidth / 4;

    // 1. Smoothing Slider
    auto s1 = advInner.removeFromLeft(itemW);
    lblSmoothing.setBounds(s1.removeFromTop(kAdvLabelH));
    s1.removeFromTop(kAdvLabelGap);
    sliderSmoothing.setBounds(s1.removeFromTop(kAdvSliderH));

    // 2. Masking Slider
    advInner.removeFromLeft(kGapW);
    auto s2 = advInner.removeFromLeft(itemW);
    lblMasking.setBounds(s2.removeFromTop(kAdvLabelH));
    s2.removeFromTop(kAdvLabelGap);
    sliderMasking.setBounds(s2.removeFromTop(kAdvSliderH));

    // 3. Whitening Slider
    advInner.removeFromLeft(kGapW);
    auto s3 = advInner.removeFromLeft(itemW);
    lblWhitening.setBounds(s3.removeFromTop(kAdvLabelH));
    s3.removeFromTop(kAdvLabelGap);
    sliderWhitening.setBounds(s3.removeFromTop(kAdvSliderH));

    // 4. Aggressiveness Slider
    advInner.removeFromLeft(kGapW);
    auto s4 = advInner;
    lblAggressiveness.setBounds(s4.removeFromTop(kAdvLabelH));
    s4.removeFromTop(kAdvLabelGap);
    sliderAggressiveness.setBounds(s4.removeFromTop(kAdvSliderH));

    area.removeFromBottom(kSectionGap);
  }

  // Main Denoising & Spectrum Canvas
  groupDenoising.setBounds(area);

  auto denoiseInner = area.reduced(kInnerPad);
  denoiseInner.removeFromTop(kGroupTitleH);

  bool isLinked = !isAdvancedVisible || btnLink.getToggleState();
  int faderBankWidth = isLinked ? kBankLinkedW : kBankUnlinkedW;

  auto faderBankArea = denoiseInner.removeFromLeft(faderBankWidth);
  lblReductionHeader.setBounds(faderBankArea.removeFromTop(kBankLabelH));
  faderBankArea.removeFromTop(kBankLabelGap);

  if (isAdvancedVisible) {
    btnLink.setBounds(faderBankArea.removeFromTop(kBankBtnH));
    faderBankArea.removeFromTop(kBankBtnGap);

    // Carve bottom area of faderBankArea for Curve split button [ Curve ][ ↺ ]
    bool isCurveOn = btnCurveToggle.getToggleState();
    auto faderBottomArea = faderBankArea.removeFromBottom(kBankBtnH);
    faderBankArea.removeFromBottom(kBankBtnGap);

    if (isCurveOn) {
      btnCurveToggle.setBounds(
          faderBottomArea.removeFromLeft(faderBottomArea.getWidth() - kCurveResetW));
      faderBottomArea.removeFromLeft(kCurveResetGap);
      btnResetCurve.setBounds(faderBottomArea);
    } else {
      btnCurveToggle.setBounds(faderBottomArea);
    }
  }

  if (isLinked) {
    sliderMasterRed.setBounds(faderBankArea);
  } else {
    int colW = (faderBankArea.getWidth() - kBankSplitGap) / 2;

    auto leftCol = faderBankArea.removeFromLeft(colW);
    faderBankArea.removeFromLeft(kBankSplitGap);
    auto rightCol = faderBankArea;

    lblMasterRed.setBounds(leftCol.removeFromTop(kBankLabelH));
    leftCol.removeFromTop(kBankLabelGap);
    sliderMasterRed.setBounds(leftCol);

    lblTonalRed.setBounds(rightCol.removeFromTop(kBankLabelH));
    rightCol.removeFromTop(kBankLabelGap);
    sliderTonalRed.setBounds(rightCol);
  }

  denoiseInner.removeFromLeft(kBankGutter);

  // Right-side Threshold (Noise Profile Offset) vertical fader bank (symmetric
  // to reduction sliders on left)
  bool isOffsetLinked = !isAdvancedVisible || btnLinkOffset.getToggleState();
  int offsetBankWidth = isOffsetLinked ? kBankLinkedW : kBankUnlinkedW;

  auto offsetBankArea = denoiseInner.removeFromRight(offsetBankWidth);
  denoiseInner.removeFromRight(kBankGutter);

  lblOffset.setBounds(offsetBankArea.removeFromTop(kBankLabelH));
  offsetBankArea.removeFromTop(kBankLabelGap);

  if (isAdvancedVisible) {
    btnLinkOffset.setBounds(offsetBankArea.removeFromTop(kBankBtnH));
    offsetBankArea.removeFromTop(kBankBtnGap);
  }

  // Carve bottom area of offsetBankArea for Advanced Controls button below
  // threshold slider
  auto offsetBottomArea = offsetBankArea.removeFromBottom(kBankBtnH);
  offsetBankArea.removeFromBottom(kBankBtnGap);
  btnAdvancedToggle.setBounds(offsetBottomArea);

  if (isOffsetLinked) {
    sliderOffset.setBounds(offsetBankArea);
  } else {
    int colW = (offsetBankArea.getWidth() - kBankSplitGap) / 2;

    auto leftCol = offsetBankArea.removeFromLeft(colW);
    offsetBankArea.removeFromLeft(kBankSplitGap);
    auto rightCol = offsetBankArea;

    lblMasterOffset.setBounds(leftCol.removeFromTop(kBankLabelH));
    leftCol.removeFromTop(kBankLabelGap);
    sliderOffset.setBounds(leftCol);

    lblTonalOffset.setBounds(rightCol.removeFromTop(kBankLabelH));
    rightCol.removeFromTop(kBankLabelGap);
    sliderTonalOffset.setBounds(rightCol);
  }

  spectralVisualizer.setBounds(denoiseInner);

  btnPreferences.toFront(false);
  footerTooltipLabel.toFront(false);
}
