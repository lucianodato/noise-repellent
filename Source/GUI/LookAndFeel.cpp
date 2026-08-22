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

#include "LookAndFeel.h"

const juce::Colour NoiseRepellentLookAndFeel::kColorNoiseProfile =
    juce::Colour(0xffe5a000);
const juce::Colour NoiseRepellentLookAndFeel::kColorDenoising =
    juce::Colour(0xff5cc0d4);
const juce::Colour NoiseRepellentLookAndFeel::kColorFineTuning =
    juce::Colour(0xffb4bfce);
const juce::Colour NoiseRepellentLookAndFeel::kColorInputSignal =
    juce::Colour(0xff6184a8);
const juce::Colour NoiseRepellentLookAndFeel::kColorTonalPeaks =
    juce::Colour(0xffe07055);
const juce::Colour NoiseRepellentLookAndFeel::kColorReductionCurve =
    juce::Colour(0xff4caf50);
const juce::Colour NoiseRepellentLookAndFeel::kColorPanelBg =
    juce::Colour(0xff343a48);
const juce::Colour NoiseRepellentLookAndFeel::kColorPanelBorder =
    juce::Colour(0xff4f586c);

NoiseRepellentLookAndFeel::NoiseRepellentLookAndFeel() {
  setColour(juce::Slider::thumbColourId, juce::Colour(0xff525c70));
  setColour(juce::Slider::trackColourId, juce::Colour(0xff1c2028));
  setColour(juce::TextButton::textColourOffId, juce::Colour(0xffd0d7e2));
  setColour(juce::TextButton::textColourOnId, kColorNoiseProfile);
  setColour(juce::ComboBox::backgroundColourId, juce::Colour(0xff262b36));
  setColour(juce::ComboBox::outlineColourId, kColorPanelBorder);
}

void NoiseRepellentLookAndFeel::drawLinearSlider(
    juce::Graphics& g, int x, int y, int width, int height, float sliderPos,
    float minSliderPos, float maxSliderPos,
    const juce::Slider::SliderStyle style, juce::Slider& slider) {
  juce::Colour fillColour =
      slider.findColour(juce::Slider::rotarySliderFillColourId, true);
  if (fillColour.isTransparent()) {
    fillColour = kColorFineTuning;
  }

  if (style == juce::Slider::LinearVertical) {
    float trackW = 4.0f;
    float trackX = x + (width - trackW) * 0.5f;

    // Background track
    g.setColour(juce::Colour(0xff1c2028));
    g.fillRoundedRectangle(trackX, (float)y, trackW, (float)height, 2.0f);

    // Solid Single-Color Fill
    g.setColour(fillColour);
    g.fillRoundedRectangle(trackX, sliderPos, trackW,
                           (float)y + (float)height - sliderPos, 2.0f);

    // Thumb Cap
    float thumbW = 18.0f;
    float thumbH = 12.0f;
    float thumbX = x + (width - thumbW) * 0.5f;
    float thumbY = sliderPos - thumbH * 0.5f;

    g.setColour(juce::Colour(0xff525c70));
    g.fillRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f);

    g.setColour(juce::Colour(0xff8492a8));
    g.drawRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f, 1.0f);

    g.setColour(fillColour);
    g.fillRect(thumbX + 5.0f, thumbY + 5.0f, 8.0f, 2.0f);
  } else if (style == juce::Slider::LinearHorizontal) {
    float trackH = 4.0f;
    float trackY = y + (height - trackH) * 0.5f;

    g.setColour(juce::Colour(0xff1c2028));
    g.fillRoundedRectangle((float)x, trackY, (float)width, trackH, 2.0f);

    g.setColour(fillColour);
    g.fillRoundedRectangle((float)x, trackY, sliderPos - (float)x, trackH,
                           2.0f);

    float thumbW = 12.0f;
    float thumbH = 16.0f;
    float thumbX = sliderPos - thumbW * 0.5f;
    float thumbY = y + (height - thumbH) * 0.5f;

    g.setColour(juce::Colour(0xff525c70));
    g.fillRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f);

    g.setColour(juce::Colour(0xff8492a8));
    g.drawRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f, 1.0f);

    g.setColour(fillColour);
    g.fillRect(thumbX + 5.0f, thumbY + 4.0f, 2.0f, 8.0f);
  }
}

void NoiseRepellentLookAndFeel::drawButtonBackground(
    juce::Graphics& g, juce::Button& button,
    const juce::Colour& backgroundColour, bool shouldDrawButtonAsHighlighted,
    bool shouldDrawButtonAsDown) {
  auto bounds = button.getLocalBounds().toFloat();
  bool isOn = button.getToggleState();

  juce::Colour base;
  if (isOn) {
    juce::Colour onCol =
        button.findColour(juce::TextButton::buttonOnColourId, false);
    base = onCol.isOpaque() ? onCol : kColorNoiseProfile;
  } else {
    juce::Colour offCol =
        button.findColour(juce::TextButton::buttonColourId, false);
    base = (offCol.isOpaque() && !offCol.isTransparent())
               ? offCol
               : (backgroundColour.isTransparent() ? juce::Colour(0xff3f4757)
                                                   : backgroundColour);
  }

  if (shouldDrawButtonAsDown)
    base = base.darker(0.2f);
  else if (shouldDrawButtonAsHighlighted)
    base = base.brighter(0.15f);

  g.setColour(base);
  g.fillRoundedRectangle(bounds, 4.0f);

  g.setColour(kColorPanelBorder);
  g.drawRoundedRectangle(bounds, 4.0f, 1.0f);
}

juce::Font NoiseRepellentLookAndFeel::getTextButtonFont(juce::TextButton&,
                                                        int) {
  return juce::FontOptions(kFontSizeLabel, juce::Font::bold);
}

void NoiseRepellentLookAndFeel::drawButtonText(
    juce::Graphics& g, juce::TextButton& button,
    bool shouldDrawButtonAsHighlighted, bool shouldDrawButtonAsDown) {
  juce::ignoreUnused(shouldDrawButtonAsHighlighted, shouldDrawButtonAsDown);
  juce::Font font(getTextButtonFont(button, button.getHeight()));
  g.setFont(font);

  bool isOn = button.getToggleState();
  juce::Colour activeBg =
      button.findColour(isOn ? juce::TextButton::buttonOnColourId
                             : juce::TextButton::buttonColourId,
                        false);

  juce::Colour textColour;
  if (activeBg == juce::Colour(0xffc0392b) ||
      activeBg == juce::Colour(0xffe74c3c)) {
    textColour = juce::Colours::white;
  } else if (activeBg == kColorNoiseProfile || isOn) {
    textColour = juce::Colour(0xff101216);
  } else {
    textColour = juce::Colour(0xffe6e8ed);
  }

  if (!button.isEnabled())
    textColour = textColour.withAlpha(0.4f);

  g.setColour(textColour);
  g.drawText(button.getButtonText(), button.getLocalBounds(),
             juce::Justification::centred);
}

void NoiseRepellentLookAndFeel::drawToggleButton(
    juce::Graphics& g, juce::ToggleButton& button,
    bool shouldDrawButtonAsHighlighted, bool shouldDrawButtonAsDown) {
  juce::ignoreUnused(shouldDrawButtonAsHighlighted, shouldDrawButtonAsDown);
  auto bounds = button.getLocalBounds().toFloat();
  bool isOn = button.getToggleState();

  g.setColour(isOn ? kColorNoiseProfile : juce::Colour(0xff3f4757));
  g.fillRoundedRectangle(bounds, 4.0f);

  g.setColour(kColorPanelBorder);
  g.drawRoundedRectangle(bounds, 4.0f, 1.0f);

  juce::Colour textColour =
      isOn ? juce::Colour(0xff101216) : juce::Colour(0xffe6e8ed);
  if (!button.isEnabled())
    textColour = textColour.withAlpha(0.4f);

  g.setColour(textColour);
  g.setFont(juce::FontOptions(kFontSizeLabel, juce::Font::bold));
  g.drawText(button.getButtonText(), bounds, juce::Justification::centred);
}

void NoiseRepellentLookAndFeel::drawComboBox(juce::Graphics& g, int width,
                                             int height, bool isButtonDown,
                                             int buttonX, int buttonY,
                                             int buttonW, int buttonH,
                                             juce::ComboBox& box) {
  auto cornerSize = 4.0f;
  juce::Rectangle<float> boxBounds(0.0f, 0.0f, (float)width, (float)height);

  g.setColour(findColour(juce::ComboBox::backgroundColourId));
  g.fillRoundedRectangle(boxBounds, cornerSize);

  g.setColour(findColour(juce::ComboBox::outlineColourId));
  g.drawRoundedRectangle(boxBounds, cornerSize, 1.0f);

  juce::Rectangle<float> arrowZone((float)buttonX, (float)buttonY,
                                   (float)buttonW, (float)buttonH);
  juce::Path path;
  path.addTriangle(arrowZone.getCentreX() - 4.0f, arrowZone.getCentreY() - 2.0f,
                   arrowZone.getCentreX() + 4.0f, arrowZone.getCentreY() - 2.0f,
                   arrowZone.getCentreX(), arrowZone.getCentreY() + 3.0f);

  g.setColour(juce::Colour(0xff808896));
  g.fillPath(path);
}

void NoiseRepellentLookAndFeel::drawGroupComponentOutline(
    juce::Graphics& g, int width, int height, const juce::String& text,
    const juce::Justification& position, juce::GroupComponent& group) {
  juce::ignoreUnused(position);

  auto textH = 14.0f;
  auto textEdgeGap = 4.0f;
  auto cs = 4.0f;

  juce::Font font(juce::FontOptions(kFontSizeLabel, juce::Font::bold));
  juce::GlyphArrangement ga;
  ga.addFittedText(font, text, 0.0f, 0.0f, 1000.0f, textH,
                   juce::Justification::left, 1);
  auto textW = text.isEmpty() ? 0.0f
                              : ga.getBoundingBox(0, -1, true).getWidth() +
                                    textEdgeGap * 2.0f;

  auto bounds =
      juce::Rectangle<float>((float)width, (float)height).reduced(0.5f);
  auto textBounds = juce::Rectangle<float>(10.0f, 0.0f, textW, textH);

  juce::Path p;
  p.startNewSubPath(textBounds.getRight(), bounds.getY() + textH * 0.5f);
  p.lineTo(bounds.getRight() - cs, bounds.getY() + textH * 0.5f);
  p.addArc(bounds.getRight() - cs * 2.0f, bounds.getY() + textH * 0.5f,
           cs * 2.0f, cs * 2.0f, 0.0f, juce::MathConstants<float>::halfPi);
  p.lineTo(bounds.getRight(), bounds.getBottom() - cs);
  p.addArc(bounds.getRight() - cs * 2.0f, bounds.getBottom() - cs * 2.0f,
           cs * 2.0f, cs * 2.0f, juce::MathConstants<float>::halfPi,
           juce::MathConstants<float>::pi);
  p.lineTo(bounds.getX() + cs, bounds.getBottom());
  p.addArc(bounds.getX(), bounds.getBottom() - cs * 2.0f, cs * 2.0f, cs * 2.0f,
           juce::MathConstants<float>::pi,
           juce::MathConstants<float>::pi * 1.5f);
  p.lineTo(bounds.getX(), bounds.getY() + textH * 0.5f + cs);
  p.addArc(bounds.getX(), bounds.getY() + textH * 0.5f, cs * 2.0f, cs * 2.0f,
           juce::MathConstants<float>::pi * 1.5f,
           juce::MathConstants<float>::twoPi);
  p.lineTo(textBounds.getX(), bounds.getY() + textH * 0.5f);

  g.setColour(group.findColour(juce::GroupComponent::outlineColourId));
  g.strokePath(p, juce::PathStrokeType(1.0f));

  if (text.isNotEmpty()) {
    g.setColour(group.findColour(juce::GroupComponent::textColourId));
    g.setFont(font);
    g.drawText(text, textBounds, juce::Justification::centred, false);
  }
}

juce::Font NoiseRepellentLookAndFeel::getComboBoxFont(juce::ComboBox&) {
  return juce::FontOptions(kFontSizeLabel, juce::Font::bold);
}

juce::Font NoiseRepellentLookAndFeel::getPopupMenuFont() {
  return juce::FontOptions(kFontSizeLabel, juce::Font::plain);
}

void NoiseRepellentLookAndFeel::drawPopupMenuSectionHeader(
    juce::Graphics& g, const juce::Rectangle<int>& area,
    const juce::String& sectionName) {
  g.setColour(juce::Colour(0xff64748b)); // Slate Gray
  g.setFont(juce::FontOptions(kFontSizeLabel - 1.0f, juce::Font::bold));
  g.drawText(sectionName, area.reduced(12, 0), juce::Justification::centredLeft,
             true);
}

void NoiseRepellentLookAndFeel::getIdealPopupMenuItemSize(
    const juce::String& text, bool isSeparator, int standardMenuItemHeight,
    int& idealWidth, int& idealHeight) {
  juce::LookAndFeel_V4::getIdealPopupMenuItemSize(
      text, isSeparator, standardMenuItemHeight, idealWidth, idealHeight);
  idealHeight = 24; // Match standard combobox height
}
