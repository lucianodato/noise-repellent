/*
noise-repellent-live -- Zero-Latency Noise Reduction JUCE Plugin

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
#include <cmath>

NoiseRepellentLiveAudioProcessorEditor::NoiseRepellentLiveAudioProcessorEditor(
    NoiseRepellentLiveAudioProcessor& p)
    : AudioProcessorEditor(&p), audioProcessor(p), spectralVisualizer(p) {
  setLookAndFeel(&customLookAndFeel);

  // ── Header ──
  addAndMakeVisible(brandLabel);
  brandLabel.setText("NOISE REPELLENT LIVE", juce::dontSendNotification);
  brandLabel.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeBrand, juce::Font::bold));
  brandLabel.setColour(juce::Label::textColourId, juce::Colour(0xffd8e0ec));
  brandLabel.setJustificationType(juce::Justification::centredLeft);

  addAndMakeVisible(btnDelta);
  btnDelta.setClickingTogglesState(true);
  btnDelta.setColour(juce::TextButton::buttonColourId,
                     juce::Colour(0xff3f4757));
  btnDelta.setColour(juce::TextButton::buttonOnColourId,
                     NoiseRepellentLookAndFeel::kColorNoiseProfile);
  btnDelta.setTooltip(
      "Overlay the reduction delta (input minus output) "
      "on the spectrum display");
  btnDelta.onClick = [this]() {
    spectralVisualizer.setDeltaVisible(btnDelta.getToggleState());
  };

  addAndMakeVisible(btnBypass);
  btnBypass.setClickingTogglesState(true);
  btnBypass.setTooltip("Engages the internal soft bypass");

  // ── Main Reduction Control ──
  addAndMakeVisible(sliderReduction);
  sliderReduction.setSliderStyle(juce::Slider::LinearVertical);
  sliderReduction.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 80, 20);
  sliderReduction.setTooltip("Amount of noise reduction applied in dB");

  addAndMakeVisible(lblReduction);
  lblReduction.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblReduction.setColour(juce::Label::textColourId, juce::Colour(0xffa8b3c4));
  lblReduction.setJustificationType(juce::Justification::centred);

  // ── Adaptive Two-Part Button (toggle + method dropdown) ──
  btnAdaptiveNoise.setClickingTogglesState(true);
  btnAdaptiveNoise.setColour(juce::TextButton::buttonColourId,
                             juce::Colour(0xff3f4757));
  btnAdaptiveNoise.setColour(juce::TextButton::buttonOnColourId,
                             NoiseRepellentLookAndFeel::kColorDenoising);
  btnAdaptiveNoise.setTooltip(
      "Continuously estimate the noise floor from the "
      "input signal. When off, the last estimate is "
      "frozen.");
  addAndMakeVisible(btnAdaptiveNoise);

  btnAdaptiveArrow.setColour(juce::TextButton::buttonColourId,
                             juce::Colour(0xff353b48));
  btnAdaptiveArrow.setTooltip("Adaptive estimation method");
  btnAdaptiveArrow.onClick = [this]() {
    juce::PopupMenu menu;
    const int currentMethod = comboMethod.getSelectedId();
    menu.addSectionHeader("ADAPTIVE ESTIMATION METHOD");
    menu.addItem(1, "SPP-MMSE (Unbiased)", true, currentMethod == 1);
    menu.addItem(2, "Brandt (Trimmed Mean)", true, currentMethod == 2);
    menu.addItem(3, "Martin (Minimum Statistics)", true, currentMethod == 3);
    menu.setLookAndFeel(&getLookAndFeel());
    menu.showMenuAsync(
        juce::PopupMenu::Options().withTargetComponent(&btnAdaptiveArrow),
        [this](int result) {
          if (result >= 1 && result <= 3) {
            comboMethod.setSelectedId(result, juce::sendNotification);
            if (auto* p =
                    audioProcessor.getAPVTS().getParameter("adaptive_method"))
              p->setValueNotifyingHost(static_cast<float>(result - 1) / 2.0f);
            if (auto* p =
                    audioProcessor.getAPVTS().getParameter("adaptive_noise"))
              p->setValueNotifyingHost(1.0f);
          }
        });
  };
  addAndMakeVisible(btnAdaptiveArrow);

  // ── Learn Button ──
  btnLearn.setClickingTogglesState(true);
  btnLearn.setColour(juce::TextButton::buttonColourId,
                     juce::Colour(0xff3f4757));
  btnLearn.setColour(juce::TextButton::buttonOnColourId,
                     NoiseRepellentLookAndFeel::kColorNoiseProfile);
  btnLearn.setTooltip(
      "Learn the noise floor shape: loop a noise-only section with "
      "this engaged and the learned profile (amber) converges on the "
      "spectrum display");
  btnLearn.onClick = [this]() {
    if (btnLearn.getToggleState()) {
      spectralVisualizer.startLearning();
    } else {
      spectralVisualizer.stopLearning();
    }
  };
  addAndMakeVisible(btnLearn);

  // ── Spectrum Display ──
  addAndMakeVisible(spectralVisualizer);

  // ── Advanced Panel ──
  addAndMakeVisible(btnAdvancedToggle);
  btnAdvancedToggle.setClickingTogglesState(true);
  btnAdvancedToggle.setColour(juce::TextButton::buttonColourId,
                              juce::Colour(0xff3f4757));
  btnAdvancedToggle.setColour(juce::TextButton::buttonOnColourId,
                              NoiseRepellentLookAndFeel::kColorNoiseProfile);
  btnAdvancedToggle.onClick = [this]() {
    isAdvancedVisible = btnAdvancedToggle.getToggleState();
    updateLayout();
  };

  addAndMakeVisible(groupAdvanced);
  groupAdvanced.setVisible(isAdvancedVisible);
  groupAdvanced.setText("ADVANCED CONTROLS");
  groupAdvanced.setColour(juce::GroupComponent::outlineColourId,
                          NoiseRepellentLookAndFeel::kColorPanelBorder);
  groupAdvanced.setColour(juce::GroupComponent::textColourId,
                          NoiseRepellentLookAndFeel::kColorFineTuning);
  groupAdvanced.setInterceptsMouseClicks(false, false);

  addAndMakeVisible(comboMethod);
  comboMethod.setVisible(false);
  comboMethod.setTooltip("Noise estimation strategy used when Adaptive is on");

  addAndMakeVisible(sliderSmoothing);
  sliderSmoothing.setVisible(isAdvancedVisible);
  sliderSmoothing.setColour(juce::Slider::rotarySliderFillColourId,
                            NoiseRepellentLookAndFeel::kColorDenoising);
  sliderSmoothing.setTooltip(
      "Temporal smoothing of the ERB band gains. "
      "Higher values adapt slower and sound steadier");

  addAndMakeVisible(sliderSuppression);
  sliderSuppression.setVisible(isAdvancedVisible);
  sliderSuppression.setColour(juce::Slider::rotarySliderFillColourId,
                              NoiseRepellentLookAndFeel::kColorDenoising);
  sliderSuppression.setTooltip("Oversubtraction aggressiveness");

  addAndMakeVisible(lblSmoothing);
  lblSmoothing.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblSmoothing.setColour(juce::Label::textColourId,
                         NoiseRepellentLookAndFeel::kColorDenoising);
  lblSmoothing.setJustificationType(juce::Justification::centred);
  lblSmoothing.setVisible(isAdvancedVisible);

  addAndMakeVisible(lblSuppression);
  lblSuppression.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblSuppression.setColour(juce::Label::textColourId,
                           NoiseRepellentLookAndFeel::kColorDenoising);
  lblSuppression.setJustificationType(juce::Justification::centred);
  lblSuppression.setVisible(isAdvancedVisible);

  // ── Footer ──
  addAndMakeVisible(footerTooltipLabel);
  footerTooltipLabel.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeTooltip, juce::Font::plain));
  footerTooltipLabel.setColour(juce::Label::textColourId,
                               juce::Colour(0xffa8b3c4));
  footerTooltipLabel.setJustificationType(juce::Justification::centredLeft);

  // Attachments
  auto& apvts = audioProcessor.getAPVTS();
  attachBypass = std::make_unique<ButtonAttachment>(apvts, "bypass", btnBypass);
  attachAdaptive = std::make_unique<ButtonAttachment>(apvts, "adaptive_noise",
                                                      btnAdaptiveNoise);
  attachReduction = std::make_unique<SliderAttachment>(
      apvts, "reduction_amount", sliderReduction);
  attachMethod = std::make_unique<ComboBoxAttachment>(apvts, "adaptive_method",
                                                      comboMethod);
  attachSmoothing = std::make_unique<SliderAttachment>(
      apvts, "smoothing_factor", sliderSmoothing);
  attachSuppression = std::make_unique<SliderAttachment>(
      apvts, "suppression_strength", sliderSuppression);

  // Footer tooltip follows hovered component
  for (auto* comp : std::array<juce::Component*, 9>{
           &sliderReduction, &btnAdaptiveNoise, &btnAdaptiveArrow, &btnLearn,
           &btnBypass, &btnDelta, &sliderSmoothing, &sliderSuppression,
           &btnAdvancedToggle}) {
    comp->addMouseListener(this, false);
  }

  setResizable(false, false);
  setResizeLimits(720, 480, 1280, 840);
  setSize(760, 520);
  updateLayout();
}

NoiseRepellentLiveAudioProcessorEditor::
    ~NoiseRepellentLiveAudioProcessorEditor() {
  setLookAndFeel(nullptr);
}

void NoiseRepellentLiveAudioProcessorEditor::mouseEnter(
    const juce::MouseEvent& event) {
  if (auto* tip =
          dynamic_cast<juce::SettableTooltipClient*>(event.eventComponent)) {
    if (tip->getTooltip().isNotEmpty()) {
      footerTooltipLabel.setText(tip->getTooltip(), juce::dontSendNotification);
    }
  }
}

void NoiseRepellentLiveAudioProcessorEditor::mouseExit(
    const juce::MouseEvent&) {
  footerTooltipLabel.setText({}, juce::dontSendNotification);
}

void NoiseRepellentLiveAudioProcessorEditor::paint(juce::Graphics& g) {
  g.fillAll(NoiseRepellentLookAndFeel::kColorPanelBg);
  g.setColour(NoiseRepellentLookAndFeel::kColorPanelBorder);
  g.drawRect(getLocalBounds(), 1.0f);
}

void NoiseRepellentLiveAudioProcessorEditor::updateLayout() {
  groupAdvanced.setVisible(isAdvancedVisible);
  sliderSmoothing.setVisible(isAdvancedVisible);
  sliderSuppression.setVisible(isAdvancedVisible);
  lblSmoothing.setVisible(isAdvancedVisible);
  lblSuppression.setVisible(isAdvancedVisible);
  resized();
}

void NoiseRepellentLiveAudioProcessorEditor::resized() {
  auto area = getLocalBounds().reduced(10);

  // ── Header ──
  auto header = area.removeFromTop(36);
  brandLabel.setBounds(header.removeFromLeft(260));
  btnBypass.setBounds(header.removeFromRight(90).reduced(0, 4));
  btnDelta.setBounds(header.removeFromRight(74).reduced(0, 4));

  area.removeFromBottom(4); // gap above footer

  // ── Footer ──
  auto footer = area.removeFromBottom(20);
  footerTooltipLabel.setBounds(footer);

  area.removeFromBottom(6); // gap above advanced panel

  // ── Advanced Panel (bottom, collapsible) ──
  const int advancedHeight = isAdvancedVisible ? 120 : 30;
  auto advancedArea = area.removeFromBottom(advancedHeight);
  btnAdvancedToggle.setBounds(
      advancedArea.removeFromTop(26).removeFromLeft(120));
  advancedArea.removeFromTop(4); // gap below toggle
  if (isAdvancedVisible) {
    groupAdvanced.setBounds(advancedArea);
    // Top inset clears the group title text
    auto inner = advancedArea.reduced(14, 8);
    inner.removeFromTop(16);
    const int half = inner.getWidth() / 2;

    auto col2 = inner.removeFromLeft(half);
    lblSmoothing.setBounds(col2.removeFromTop(16));
    const int knob2 = std::min(col2.getWidth() - 40, col2.getHeight() - 4);
    sliderSmoothing.setBounds(col2.withSizeKeepingCentre(knob2, knob2));

    auto col3 = inner;
    lblSuppression.setBounds(col3.removeFromTop(16));
    const int knob3 = std::min(col3.getWidth() - 40, col3.getHeight() - 4);
    sliderSuppression.setBounds(col3.withSizeKeepingCentre(knob3, knob3));
  }

  area.removeFromBottom(6); // gap above main area

  // ── Main Area: Reduction Slider (left) + Spectrum (rest) ──
  auto main = area;
  auto reductionCol = main.removeFromLeft(112);
  main.removeFromLeft(8); // gap between column and spectrum
  lblReduction.setBounds(reductionCol.removeFromTop(18));
  reductionCol.removeFromTop(6);

  // Learn button (full width), then Adaptive + dropdown arrow row
  auto learnBounds = reductionCol.removeFromTop(28).reduced(12, 2);
  btnLearn.setBounds(learnBounds);
  reductionCol.removeFromTop(4);
  auto adaptRow = reductionCol.removeFromTop(28).reduced(12, 2);
  btnAdaptiveNoise.setBounds(adaptRow.removeFromLeft(adaptRow.getWidth() - 22));
  adaptRow.removeFromLeft(4);
  btnAdaptiveArrow.setBounds(adaptRow.removeFromRight(18));

  reductionCol.removeFromTop(8);
  sliderReduction.setBounds(reductionCol);

  spectralVisualizer.setBounds(main);
}
