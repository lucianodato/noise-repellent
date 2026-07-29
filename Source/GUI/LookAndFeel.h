/*
noise-repellent -- Noise Reduction JUCE Plugin

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.
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

    void drawToggleButton(juce::Graphics& g, juce::ToggleButton& button,
                          bool shouldDrawButtonAsHighlighted,
                          bool shouldDrawButtonAsDown) override;

    // Domain Colors
    static const juce::Colour kColorNoiseProfile; // Warm Amber 0xffe5a000
    static const juce::Colour kColorDenoising;    // Slate Blue  0xff00b4d8
    static const juce::Colour kColorFineTuning;   // Slate Gray  0xff64748b
    static const juce::Colour kColorInputSignal;  // Soft Teal   0xff486581
    static const juce::Colour kColorTonalPeaks;   // Muted Violet 0xff9d65c9
    static const juce::Colour kColorPanelBg;      // Matte Dark  0xff171a21
    static const juce::Colour kColorPanelBorder;  // Border      0xff272c37
};
