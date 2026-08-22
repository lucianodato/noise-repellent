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

NoiseRepellentAudioProcessorEditor::NoiseRepellentAudioProcessorEditor(
    NoiseRepellentAudioProcessor& p)
    : AudioProcessorEditor(&p), audioProcessor(p), spectralVisualizer(p) {
  setLookAndFeel(&customLookAndFeel);

  // 1. Native Resizable Window Settings
  setResizable(true, true);
  setResizeLimits(840, 480, 1600, 1000);
  setSize(940, 560);

  // Brand Header
  brandLabel.setText("NOISE REPELLENT", juce::dontSendNotification);
  brandLabel.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeBrand, juce::Font::bold));
  brandLabel.setColour(juce::Label::textColourId, juce::Colours::white);
  addAndMakeVisible(brandLabel);

  // Preferences Header Button
  addAndMakeVisible(btnPreferences);
  btnPreferences.setTooltip("Preferences");
  btnPreferences.setColour(juce::TextButton::buttonColourId,
                           juce::Colour(0xff3a4150));
  btnPreferences.onClick = [this]() {
    auto* tooltipParam =
        audioProcessor.getAPVTS().getParameter("show_tooltips");
    auto* hudParam = audioProcessor.getAPVTS().getParameter("show_hud");

    bool tooltipVal =
        (tooltipParam != nullptr && tooltipParam->getValue() > 0.5f);
    bool hudVal = (hudParam != nullptr && hudParam->getValue() > 0.5f);

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
      {"1D Spectral (Fast & Low CPU)", "2D NLM Patch (High Quality)"}, 1);
  addAndMakeVisible(comboAlgoMode);
  comboAlgoMode.onChange = [this]() { updateLayout(); };

  // Header Action Toggles
  isAdvancedVisible =
      audioProcessor.getAPVTS().getRawParameterValue("show_advanced")->load() >
      0.5f;
  btnAdvancedToggle.setClickingTogglesState(true);
  addAndMakeVisible(btnAdvancedToggle);
  btnAdvancedToggle.setColour(juce::TextButton::buttonColourId,
                              juce::Colour(0xff3f4757));
  btnAdvancedToggle.setColour(juce::TextButton::buttonOnColourId,
                              NoiseRepellentLookAndFeel::kColorNoiseProfile);
  btnAdvancedToggle.onClick = [this]() {
    isAdvancedVisible = btnAdvancedToggle.getToggleState();
    updateLayout();
  };

  btnDelta.setClickingTogglesState(true);
  btnDelta.setColour(juce::TextButton::buttonColourId,
                     juce::Colour(0xff3f4757));
  btnDelta.setColour(juce::TextButton::buttonOnColourId,
                     NoiseRepellentLookAndFeel::kColorNoiseProfile);
  addAndMakeVisible(btnDelta);
  btnDelta.onClick = [this]() { updateLayout(); };

  btnBypass.setClickingTogglesState(true);
  btnBypass.setColour(juce::TextButton::buttonColourId,
                      juce::Colour(0xff3f4757));
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
      juce::Colour(0xffc0392b)); // Prominent Studio Crimson Red CTA
  btnLearn.setColour(juce::TextButton::buttonOnColourId,
                     juce::Colour(0xffe74c3c)); // Active Learning Red
  btnLearn.onClick = [this]() { updateProfileStatus(); };

  btnAdaptiveNoise.setClickingTogglesState(true);
  btnAdaptiveNoise.setButtonText("Adaptive");
  btnAdaptiveNoise.setTooltip(
      "Toggle continuous adaptive noise estimation.\nWorks standalone or on "
      "top of a manually learned profile.");
  addAndMakeVisible(btnAdaptiveNoise);
  btnAdaptiveNoise.setColour(juce::TextButton::buttonColourId,
                             juce::Colour(0xff3f4757));
  btnAdaptiveNoise.setColour(juce::TextButton::buttonOnColourId,
                             NoiseRepellentLookAndFeel::kColorDenoising);
  btnAdaptiveNoise.onClick = [this]() {
    updateLayout();
    updateProfileStatus();
  };

  addAndMakeVisible(btnAdaptiveArrow);
  const juce::String adaptiveMethodTip =
      "Select Adaptive Estimation Method: SPP-MMSE (best for speech & dynamic "
      "noise),\nBrandt (best for steady hiss & fans), or Martin (best for slow "
      "background).";
  btnAdaptiveArrow.setTooltip(adaptiveMethodTip);
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
  btnResetProfile.setTooltip("Reset noise profile");
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
  btnLink.setButtonText(btnLink.getToggleState() ? "Linked" : "Unlinked");
  btnLink.onClick = [this]() {
    bool isLinked = btnLink.getToggleState();
    btnLink.setButtonText(isLinked ? "Linked" : "Unlinked");
    if (isLinked) {
      // Reset reduction parameters to default (15.0 dB) when relinked
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
  btnCurveToggle.setTooltip("Enable per-frequency reduction bias curve.");
  addAndMakeVisible(btnCurveToggle);
  btnCurveToggle.onClick = [this]() { updateLayout(); };

  btnResetCurve.setTooltip("Reset reduction curve to flat 0 dB line.");
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
                         NoiseRepellentLookAndFeel::kColorDenoising);
  addAndMakeVisible(sliderOffset);
  addAndMakeVisible(lblOffset);

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
  sliderSuppression.setColour(juce::Slider::rotarySliderFillColourId,
                              NoiseRepellentLookAndFeel::kColorDenoising);

  addAndMakeVisible(sliderSmoothing);
  addAndMakeVisible(sliderMasking);
  addAndMakeVisible(sliderWhitening);
  addAndMakeVisible(sliderSuppression);

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

  lblSuppression.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblSuppression.setColour(juce::Label::textColourId,
                           NoiseRepellentLookAndFeel::kColorDenoising);
  lblSuppression.setJustificationType(juce::Justification::centred);
  addAndMakeVisible(lblSuppression);

  // Advanced Controls Toggle Button
  btnAdvancedToggle.setClickingTogglesState(true);
  btnAdvancedToggle.setButtonText("ADVANCED");
  btnAdvancedToggle.setTooltip(
      "Toggle Advanced DSP Controls (Smoothing, Masking, Whitening, "
      "Suppression, Profile Morphing).");
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
  attachLink =
      std::make_unique<ButtonAttachment>(apvts, "link_reduction", btnLink);
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
  attachAggressiveness = std::make_unique<SliderAttachment>(
      apvts, "aggressiveness", sliderAggressiveness);
  attachSmoothing = std::make_unique<SliderAttachment>(
      apvts, "smoothing_factor", sliderSmoothing);
  attachMasking =
      std::make_unique<SliderAttachment>(apvts, "masking_depth", sliderMasking);
  attachWhitening = std::make_unique<SliderAttachment>(
      apvts, "whitening_factor", sliderWhitening);
  attachSuppression = std::make_unique<SliderAttachment>(
      apvts, "suppression_strength", sliderSuppression);

  // Control Tooltip Descriptions
  btnPreferences.setTooltip("Plugin preferences menu.");
  comboAlgoMode.setTooltip(
      "Select denoising algorithm: 1D STFT (fast & low latency)\nor 2D NLM "
      "Patch (high quality non-local means processing).");
  lblAlgoHeader.setTooltip(
      "Select denoising algorithm: 1D STFT (fast & low latency)\nor 2D NLM "
      "Patch (high quality non-local means processing).");
  btnAdvancedToggle.setTooltip(
      "Show Advanced DSP Controls (Smoothing, Masking Protect, Whitening,\n "
      "Suppression, Bias Curve, Tonal Split & Aggressiveness).");
  btnLearn.setTooltip(
      "Capture static noise profile from current audio input\n(supports "
      "multi-section accumulation).");
  btnResetProfile.setTooltip("Reset noise profile and clear learned data.");
  btnAdaptiveNoise.setTooltip(
      "Toggle continuous adaptive noise estimation.\nWorks standalone or on "
      "top of a manually learned profile.");
  btnAdaptiveArrow.setTooltip(adaptiveMethodTip);
  comboMethod.setTooltip(adaptiveMethodTip);
  lblMethod.setTooltip(adaptiveMethodTip);

  const juce::String masterRedTip =
      "Adjust noise reduction level in decibels\n(0 to 40 dB across all "
      "bands).";
  const juce::String tonalRedTip =
      "Adjust reduction level for tonal noise components\n(0 to 40 dB for "
      "harmonic peaks).";
  sliderMasterRed.setTooltip(masterRedTip);
  lblMasterRed.setTooltip(masterRedTip);
  lblReductionHeader.setTooltip(masterRedTip);
  sliderTonalRed.setTooltip(tonalRedTip);
  lblTonalRed.setTooltip(tonalRedTip);

  sliderOffset.setTooltip(
      "Shift noise profile threshold up or down\nin decibels (-12 to +12 dB).");
  lblOffset.setTooltip(
      "Shift noise profile threshold up or down\nin decibels (-12 to +12 dB).");
  const juce::String aggrTip =
      "Adjust noise profile aggressiveness using statistical variance.\n"
      "(-1.0: Median, 0.0: Mean, +1.0: Mean + 2 Standard Deviations).";
  sliderAggressiveness.setTooltip(aggrTip);
  lblAggressiveness.setTooltip(aggrTip);
  btnDelta.setTooltip(
      "Listen to removed noise signal\n(residual audio output).");
  btnBypass.setTooltip(
      "Soft bypass plugin processing\nwith smooth wet/dry fade transition.");

  sliderSmoothing.setTooltip(
      "Apply temporal smoothing across spectral frames\nto reduce musical "
      "noise bubbling artifacts.");
  sliderMasking.setTooltip(
      "Adjust psychoacoustic masking threshold\nto protect quiet musical "
      "transients.");
  sliderWhitening.setTooltip(
      "Spectral whitening factor to equalize\nthe residual noise floor "
      "spectrum.");
  sliderSuppression.setTooltip(
      "Maximum noise suppression floor (-dB)\nto prevent over-attenuation "
      "artifacts.");
  btnLink.setTooltip(
      "Link broadband and tonal reduction controls together\nfor unified "
      "adjustment.");
  btnCurveToggle.setTooltip(
      "Enable per-frequency reduction bias curve\noverlay on spectral "
      "display.");
  btnResetCurve.setTooltip("Reset reduction curve to flat 0 dB line.");

  // Initial Layout update
  isAdvancedVisible = btnAdvancedToggle.getToggleState();
  updateLayout();
  updateSliderLabels();
  startTimerHz(30);

  // Initial default instruction text
  footerTooltipLabel.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeTooltip, juce::Font::plain));
  footerTooltipLabel.setColour(juce::Label::textColourId,
                               juce::Colour(0xff94a3b8));
  footerTooltipLabel.setJustificationType(juce::Justification::left);
  addAndMakeVisible(footerTooltipLabel);

  auto* showTooltipsParam =
      audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
  if (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f) {
    footerTooltipLabel.setText(
        "Loop a noise-only segment in your DAW and click Learn Noise\nto "
        "capture a profile, then adjust reduction level.",
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
      "Loop a noise-only segment in your DAW and click Learn Noise\nto capture "
      "a profile, then adjust reduction level.",
      juce::dontSendNotification);
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
        "Loop a noise-only segment in your DAW and click Learn Noise\nto "
        "capture a profile, then adjust reduction level.",
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
        "Loop a noise-only segment in your DAW and click Learn Noise\nto "
        "capture a profile, then adjust reduction level.",
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
  sliderSuppression.setVisible(isAdvancedVisible);
  lblSuppression.setVisible(isAdvancedVisible);

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
  lblSuppression.setText(
      "SUPPRESSION: " +
          juce::String(static_cast<int>(sliderSuppression.getValue())) + "%",
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
                       juce::Colour(0xffe74c3c)); // Active Learning Red
    btnLearn.setColour(juce::TextButton::buttonOnColourId,
                       juce::Colour(0xffe74c3c));
    btnLearn.setTooltip(
        "Capturing noise profile from current audio playback...");
  } else if (hasProfile) {
    btnLearn.setButtonText("+ Learn");
    btnLearn.setColour(juce::TextButton::buttonColourId,
                       NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btnLearn.setColour(juce::TextButton::buttonOnColourId,
                       NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btnLearn.setTooltip(
        "Capture and accumulate an additional noise section\ninto the current "
        "profile.");
  } else {
    btnLearn.setButtonText("Learn Noise");
    btnLearn.setColour(juce::TextButton::buttonColourId,
                       juce::Colour(0xffc0392b)); // Prominent Crimson Red CTA
    btnLearn.setColour(juce::TextButton::buttonOnColourId,
                       juce::Colour(0xffe74c3c)); // Active Learning Red
    btnLearn.setTooltip(
        "Loop a noise-only segment in your DAW and click\nto capture a noise "
        "profile.");
  }
  btnLearn.repaint();

  if (isAdaptive) {
    btnAdaptiveNoise.setColour(juce::TextButton::buttonColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
    btnAdaptiveNoise.setColour(juce::TextButton::buttonOnColourId,
                               NoiseRepellentLookAndFeel::kColorDenoising);
  } else {
    btnAdaptiveNoise.setColour(juce::TextButton::buttonColourId,
                               juce::Colour(0xff3f4757));
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
  sliderOffset.setTooltip(
      "Shift noise profile threshold up or down\nin decibels (-12 to +12 dB).");
  lblOffset.setTooltip(
      "Shift noise profile threshold up or down\nin decibels (-12 to +12 dB).");

  // Aggressiveness (Profile Morphing) enabled in Noise Profile box when a
  // manual profile exists and Advanced Controls is ON
  bool aggressivenessEnabled = pluginActive && isAdvancedVisible && hasProfile;
  sliderAggressiveness.setEnabled(aggressivenessEnabled);
  lblAggressiveness.setEnabled(aggressivenessEnabled);

  if (aggressivenessEnabled) {
    const juce::String aggrTip =
        "Adjust noise profile aggressiveness using statistical variance.\n"
        "(-1.0: Median, 0.0: Mean, +1.0: Mean + 2 Standard Deviations).";
    sliderAggressiveness.setTooltip(aggrTip);
    lblAggressiveness.setTooltip(aggrTip);
  } else {
    sliderAggressiveness.setTooltip(
        "Profile Morphing is available in Advanced Controls\nwhen a manual "
        "learned noise profile is captured.");
    lblAggressiveness.setTooltip(
        "Profile Morphing is available in Advanced Controls\nwhen a manual "
        "learned noise profile is captured.");
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
    btnLink.setTooltip(
        "Link broadband and tonal reduction controls together\nfor unified "
        "adjustment.");
  } else {
    btnLink.setTooltip(
        "Unlinking tonal reduction is disabled in Standalone Adaptive "
        "mode\n(requires a captured manual profile to detect tonal peaks).");
  }

  lblReductionHeader.setEnabled(pluginActive);
  sliderMasterRed.setEnabled(pluginActive);
  lblMasterRed.setEnabled(pluginActive);

  bool isLinked = btnLink.getToggleState();
  bool tonalEnabled = pluginActive && !isLinked && allowUnlink;
  sliderTonalRed.setEnabled(tonalEnabled);
  lblTonalRed.setEnabled(tonalEnabled);

  const juce::String masterRedTip =
      isLinked
          ? "Adjust noise reduction level in decibels\n(0 to 40 dB across all "
            "bands)."
          : "Adjust broadband noise reduction level in decibels\n(0 to 40 dB).";
  const juce::String tonalRedTip =
      "Adjust reduction level for tonal noise components\n(0 to 40 dB).";

  lblReductionHeader.setTooltip(masterRedTip);
  sliderMasterRed.setTooltip(masterRedTip);
  lblMasterRed.setTooltip(masterRedTip);
  sliderTonalRed.setTooltip(tonalRedTip);
  lblTonalRed.setTooltip(tonalRedTip);

  btnCurveToggle.setEnabled(pluginActive);
  btnResetCurve.setEnabled(pluginActive && btnCurveToggle.getToggleState());

  // Advanced controls
  sliderSmoothing.setEnabled(pluginActive);
  sliderMasking.setEnabled(pluginActive);
  sliderWhitening.setEnabled(pluginActive);
  sliderSuppression.setEnabled(pluginActive);
  lblSmoothing.setEnabled(pluginActive);
  lblMasking.setEnabled(pluginActive);
  lblWhitening.setEnabled(pluginActive);
  lblSuppression.setEnabled(pluginActive);
  comboMethod.setEnabled(pluginActive);
  lblMethod.setEnabled(pluginActive);
  groupAdvanced.setEnabled(pluginActive);

  bool is2D = (comboAlgoMode.getSelectedItemIndex() == 1);
  const juce::String smoothingTip =
      is2D ? "Adjust NLM patch similarity filtering strength (h-parameter)\nto "
             "control 2D time-frequency smoothing."
           : "Apply temporal smoothing across spectral frames\nto reduce "
             "musical noise bubbling artifacts.";

  sliderSmoothing.setTooltip(smoothingTip);
  lblSmoothing.setTooltip(smoothingTip);

  // Status label (HUD overlay on FFT spectrum chart)
  bool isDelta = btnDelta.getToggleState();

  if (isBypassed) {
    lblProfileStatus.setText("STATUS: BYPASSED", juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               juce::Colour(0xff808896));
  } else if (isDelta) {
    lblProfileStatus.setText("STATUS: DELTA (RESIDUAL LISTEN)",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               NoiseRepellentLookAndFeel::kColorTonalPeaks);
  } else if (isLearning) {
    lblProfileStatus.setText("STATUS: CAPTURING NOISE PROFILE...",
                             juce::dontSendNotification);
    lblProfileStatus.setColour(juce::Label::textColourId,
                               juce::Colour(0xffe74c3c));
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
                               juce::Colour(0xff808896));
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
    g.setColour(juce::Colour(0xff4f586c));

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
    if (sliderWhitening.getWidth() > 0 && sliderSuppression.getWidth() > 0) {
      float x2 =
          (float)(sliderWhitening.getRight() + sliderSuppression.getX()) / 2.0f;
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
  g.setColour(juce::Colour(0xff94a3b8));
  g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                              juce::Font::plain));
  g.drawText(
      "DSP processing in progress...\nVisualizer and UI disabled for maximum "
      "speed.",
      contentArea, juce::Justification::centred, false);
}

void NoiseRepellentAudioProcessorEditor::resized() {
  auto area = getLocalBounds().reduced(12);

  // Header (Fixed 54px with spacious headroom for profile controls)
  auto headerArea = area.removeFromTop(54);

  // Left Title & Options Block
  brandLabel.setBounds(
      headerArea.removeFromLeft(135).withSizeKeepingCentre(135, 26));
  headerArea.removeFromLeft(4);
  btnPreferences.setBounds(
      headerArea.removeFromLeft(22).withSizeKeepingCentre(22, 22));

  // Right Action Buttons Block
  btnBypass.setBounds(
      headerArea.removeFromRight(60).withSizeKeepingCentre(60, 24));
  headerArea.removeFromRight(6);
  btnDelta.setBounds(
      headerArea.removeFromRight(50).withSizeKeepingCentre(50, 24));

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
  int kProfileWidth =
      isAdvancedVisible ? (buttonsWidth + 12 + 130) : (buttonsWidth + 12);
  int combinedHeaderWidth = kAlgoWidth + kHeaderGap + kProfileWidth;

  if (availMiddleWidth > combinedHeaderWidth) {
    int leftPad = (availMiddleWidth - combinedHeaderWidth) / 2;
    headerArea.removeFromLeft(leftPad);
  }

  // Processing Engine Column in Header
  auto algoCol = headerArea.removeFromLeft(kAlgoWidth);
  int algoYPad = std::max(0, (algoCol.getHeight() - 44) / 2);
  algoCol.removeFromTop(algoYPad);
  lblAlgoHeader.setBounds(algoCol.removeFromTop(15));
  algoCol.removeFromTop(3);
  comboAlgoMode.setBounds(algoCol.removeFromTop(26));

  headerArea.removeFromLeft(kHeaderGap);

  // Encapsulated Noise Profile Group Box in Header
  auto profileGroupArea = headerArea.removeFromLeft(kProfileWidth);
  groupProfile.setBounds(profileGroupArea);

  auto profileInner = profileGroupArea.reduced(6, 4);
  profileInner.removeFromTop(10);

  auto profileRow = profileInner;

  auto bLearn = profileRow.removeFromLeft(kLearnW);
  btnLearn.setBounds(bLearn.withSizeKeepingCentre(kLearnW, 24));
  profileRow.removeFromLeft(kBtnGap);

  auto bAdapt = profileRow.removeFromLeft(kAdaptW);
  btnAdaptiveNoise.setBounds(bAdapt.withSizeKeepingCentre(kAdaptW, 24));

  if (isAdvancedVisible) {
    auto bAdaptArrow = profileRow.removeFromLeft(kArrowW);
    btnAdaptiveArrow.setBounds(bAdaptArrow.withSizeKeepingCentre(kArrowW, 24));
  }
  profileRow.removeFromLeft(kBtnGap);

  auto bReset = profileRow.removeFromLeft(kResetW);
  btnResetProfile.setBounds(bReset.withSizeKeepingCentre(kResetW, 24));
  profileRow.removeFromLeft(8);

  if (isAdvancedVisible) {
    // Aggressiveness Slider Stack (Vertically Centered on right side of Noise
    // Profile box)
    auto aggrCol = profileRow;
    constexpr int kStackHeight = 14 + 3 + 16;
    int aggrYPad = std::max(0, (aggrCol.getHeight() - kStackHeight) / 2);
    aggrCol.removeFromTop(aggrYPad);

    lblAggressiveness.setBounds(aggrCol.removeFromTop(14));
    aggrCol.removeFromTop(3);
    sliderAggressiveness.setBounds(aggrCol.removeFromTop(16));
  }

  area.removeFromTop(8);

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
    area.removeFromBottom(8);

    if (showTooltips && showHud) {
      auto leftFooter =
          footerArea.removeFromLeft(footerArea.getWidth() * 55 / 100);
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
  // Whitening, Suppression)
  if (isAdvancedVisible) {
    auto advArea = area.removeFromBottom(68);
    groupAdvanced.setBounds(advArea);

    auto advInner = advArea.reduced(10, 4);
    advInner.removeFromTop(12);

    constexpr int kNumGaps = 3;
    constexpr int kGapW = 14;
    int totalAvailWidth = advInner.getWidth() - (kNumGaps * kGapW);
    int itemW = totalAvailWidth / 4;

    // 1. Smoothing Slider
    auto s1 = advInner.removeFromLeft(itemW);
    lblSmoothing.setBounds(s1.removeFromTop(16));
    s1.removeFromTop(2);
    sliderSmoothing.setBounds(s1.removeFromTop(20));

    // 2. Masking Slider
    advInner.removeFromLeft(kGapW);
    auto s2 = advInner.removeFromLeft(itemW);
    lblMasking.setBounds(s2.removeFromTop(16));
    s2.removeFromTop(2);
    sliderMasking.setBounds(s2.removeFromTop(20));

    // 3. Whitening Slider
    advInner.removeFromLeft(kGapW);
    auto s3 = advInner.removeFromLeft(itemW);
    lblWhitening.setBounds(s3.removeFromTop(16));
    s3.removeFromTop(2);
    sliderWhitening.setBounds(s3.removeFromTop(20));

    // 4. Suppression Slider
    advInner.removeFromLeft(kGapW);
    auto s4 = advInner;
    lblSuppression.setBounds(s4.removeFromTop(16));
    s4.removeFromTop(2);
    sliderSuppression.setBounds(s4.removeFromTop(20));

    area.removeFromBottom(8);
  }

  // Main Denoising & Spectrum Canvas
  groupDenoising.setBounds(area);

  auto denoiseInner = area.reduced(10);
  denoiseInner.removeFromTop(10);

  bool isLinked = !isAdvancedVisible || btnLink.getToggleState();
  int faderBankWidth = isLinked ? 95 : 155;

  auto faderBankArea = denoiseInner.removeFromLeft(faderBankWidth);
  lblReductionHeader.setBounds(faderBankArea.removeFromTop(16));
  faderBankArea.removeFromTop(2);

  if (isAdvancedVisible) {
    btnLink.setBounds(faderBankArea.removeFromTop(22));
    faderBankArea.removeFromTop(6);

    // Carve bottom area of faderBankArea for Curve split button [ Curve ][ ↺ ]
    bool isCurveOn = btnCurveToggle.getToggleState();
    auto faderBottomArea = faderBankArea.removeFromBottom(22);
    faderBankArea.removeFromBottom(6);

    if (isCurveOn) {
      btnCurveToggle.setBounds(
          faderBottomArea.removeFromLeft(faderBottomArea.getWidth() - 22));
      faderBottomArea.removeFromLeft(2);
      btnResetCurve.setBounds(faderBottomArea);
    } else {
      btnCurveToggle.setBounds(faderBottomArea);
    }
  }

  if (isLinked) {
    sliderMasterRed.setBounds(faderBankArea);
  } else {
    int colW = (faderBankArea.getWidth() - 6) / 2;

    auto leftCol = faderBankArea.removeFromLeft(colW);
    faderBankArea.removeFromLeft(6);
    auto rightCol = faderBankArea;

    lblMasterRed.setBounds(leftCol.removeFromTop(16));
    leftCol.removeFromTop(2);
    sliderMasterRed.setBounds(leftCol);

    lblTonalRed.setBounds(rightCol.removeFromTop(16));
    rightCol.removeFromTop(2);
    sliderTonalRed.setBounds(rightCol);
  }

  denoiseInner.removeFromLeft(10);

  // Right-side Threshold (Noise Profile Offset) vertical fader bank (symmetric
  // to reduction sliders on left)
  auto offsetBankArea = denoiseInner.removeFromRight(95);
  denoiseInner.removeFromRight(10);

  lblOffset.setBounds(offsetBankArea.removeFromTop(16));
  offsetBankArea.removeFromTop(2);

  // Carve bottom area of offsetBankArea for Advanced Controls button below
  // threshold slider
  auto offsetBottomArea = offsetBankArea.removeFromBottom(22);
  offsetBankArea.removeFromBottom(6);
  btnAdvancedToggle.setBounds(offsetBottomArea);

  sliderOffset.setBounds(offsetBankArea);

  spectralVisualizer.setBounds(denoiseInner);

  btnPreferences.toFront(false);
  footerTooltipLabel.toFront(false);
}
