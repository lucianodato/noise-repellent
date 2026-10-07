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
      "Monitor the removed noise: outputs dry minus denoised "
      "instead of the denoised signal");
  btnDelta.onClick = [this]() {
    audioProcessor.setDeltaMonitoring(btnDelta.getToggleState());
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

  // ── Learn Toggle (sticky: tracker converges while engaged) ──
  btnLearn.setClickingTogglesState(true);
  btnLearn.setColour(juce::TextButton::buttonColourId,
                     juce::Colour(0xff3f4757));
  btnLearn.setColour(juce::TextButton::buttonOnColourId,
                     NoiseRepellentLookAndFeel::kColorNoiseProfile);
  btnLearn.setTooltip(
      "Learn the noise profile: loop a noise-only section while "
      "engaged. Disengage to freeze the captured threshold (amber)");
  addAndMakeVisible(btnLearn);

  // ── Spectrum Display ──
  addAndMakeVisible(spectralVisualizer);

  // ── Threshold (smoothing is automatic at tuned defaults) ──
  addAndMakeVisible(sliderThreshold);
  sliderThreshold.setSliderStyle(juce::Slider::LinearVertical);
  sliderThreshold.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 80, 20);
  sliderThreshold.setTooltip(
      "Gate threshold offset in dB, like the full denoiser. "
      "Higher removes more noise, lower passes more through");

  addAndMakeVisible(lblThreshold);
  lblThreshold.setFont(juce::FontOptions(
      NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
  lblThreshold.setColour(juce::Label::textColourId,
                         NoiseRepellentLookAndFeel::kColorDenoising);
  lblThreshold.setJustificationType(juce::Justification::centred);

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
  attachLearn = std::make_unique<ButtonAttachment>(apvts, "learning", btnLearn);
  attachReduction = std::make_unique<SliderAttachment>(
      apvts, "reduction_amount", sliderReduction);
  attachThreshold = std::make_unique<SliderAttachment>(apvts, "threshold_db",
                                                       sliderThreshold);

  // Footer tooltip follows hovered component
  for (auto* comp : std::array<juce::Component*, 5>{
           &sliderReduction, &btnLearn, &btnBypass, &btnDelta,
           &sliderThreshold}) {
    comp->addMouseListener(this, false);
  }

  setResizable(false, false);
  setResizeLimits(720, 480, 1280, 840);
  setSize(760, 560);
  resized();
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

  area.removeFromBottom(6); // gap above main area

  // ── Main Area: Reduction (left) + Spectrum (center) + Threshold
  // (right, like the full denoiser offset bank) ──
  auto main = area;
  auto reductionCol = main.removeFromLeft(112);
  main.removeFromLeft(8); // gap between column and spectrum
  auto thresholdCol = main.removeFromRight(95);
  main.removeFromRight(8); // gap between spectrum and threshold
  lblReduction.setBounds(reductionCol.removeFromTop(18));
  reductionCol.removeFromTop(6);

  // Learn toggle (full width)
  btnLearn.setBounds(reductionCol.removeFromTop(28).reduced(12, 2));

  reductionCol.removeFromTop(8);
  sliderReduction.setBounds(reductionCol);

  lblThreshold.setBounds(thresholdCol.removeFromTop(18));
  thresholdCol.removeFromTop(6);
  sliderThreshold.setBounds(thresholdCol);

  spectralVisualizer.setBounds(main);
}
