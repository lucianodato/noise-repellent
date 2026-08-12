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

#include <juce_gui_basics/juce_gui_basics.h>

class NoiseRepellentLookAndFeel : public juce::LookAndFeel_V4
{
public:
    NoiseRepellentLookAndFeel();
    ~NoiseRepellentLookAndFeel() override = default;

    void drawLinearSlider(juce::Graphics& g, int x, int y, int width, int height,
                          float sliderPos, float minSliderPos, float maxSliderPos,
                          const juce::Slider::SliderStyle style, juce::Slider& slider) override;

    void drawButtonBackground(juce::Graphics& g, juce::Button& button,
                               const juce::Colour& backgroundColour,
                               bool shouldDrawButtonAsHighlighted,
                               bool shouldDrawButtonAsDown) override;

    juce::Font getTextButtonFont(juce::TextButton&, int buttonHeight) override;

    void drawButtonText(juce::Graphics& g, juce::TextButton& button,
                        bool shouldDrawButtonAsHighlighted,
                        bool shouldDrawButtonAsDown) override;

    void drawToggleButton(juce::Graphics& g, juce::ToggleButton& button,
                          bool shouldDrawButtonAsHighlighted,
                          bool shouldDrawButtonAsDown) override;

    void drawComboBox(juce::Graphics& g, int width, int height, bool isButtonDown,
                      int buttonX, int buttonY, int buttonW, int buttonH,
                      juce::ComboBox& box) override;

    void drawGroupComponentOutline(juce::Graphics& g, int width, int height,
                                   const juce::String& text, const juce::Justification& position,
                                   juce::GroupComponent& group) override;

    juce::Font getComboBoxFont(juce::ComboBox& box) override;
    juce::Font getPopupMenuFont() override;
    void drawPopupMenuSectionHeader(juce::Graphics& g, const juce::Rectangle<int>& area, const juce::String& sectionName) override;
    void getIdealPopupMenuItemSize(const juce::String& text, bool isSeparator, int standardMenuItemHeight, int& idealWidth, int& idealHeight) override;

    // Domain Colors
    static const juce::Colour kColorNoiseProfile; // Warm Amber 0xffe5a000
    static const juce::Colour kColorDenoising;    // Slate Blue  0xff00b4d8
    static const juce::Colour kColorFineTuning;   // Slate Gray  0xff64748b
    static const juce::Colour kColorInputSignal;  // Soft Teal   0xff486581
    static const juce::Colour kColorTonalPeaks;   // Muted Violet 0xff9d65c9
    static const juce::Colour kColorReductionCurve; // Soft Green 0xff4caf50
    static const juce::Colour kColorPanelBg;      // Matte Dark  0xff171a21
    static const juce::Colour kColorPanelBorder;  // Border      0xff272c37

    // Typography Scale
    static constexpr float kFontSizeLabel = 12.0f;  // Standardized UI font size for all labels, headers, buttons, dropdowns
    static constexpr float kFontSizeTooltip = 13.5f; // Increased font size for footer tooltips
    static constexpr float kFontSizeBrand = 18.0f;  // Brand name only
};
