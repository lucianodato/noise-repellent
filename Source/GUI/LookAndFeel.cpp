/*
noise-repellent -- Noise Reduction JUCE Plugin

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.
*/

#include "LookAndFeel.h"

const juce::Colour NoiseRepellentLookAndFeel::kColorNoiseProfile = juce::Colour(0xffe5a000);
const juce::Colour NoiseRepellentLookAndFeel::kColorDenoising    = juce::Colour(0xff00b4d8);
const juce::Colour NoiseRepellentLookAndFeel::kColorFineTuning   = juce::Colour(0xff64748b);
const juce::Colour NoiseRepellentLookAndFeel::kColorInputSignal  = juce::Colour(0xff486581);
const juce::Colour NoiseRepellentLookAndFeel::kColorTonalPeaks   = juce::Colour(0xff9d65c9);
const juce::Colour NoiseRepellentLookAndFeel::kColorPanelBg      = juce::Colour(0xff171a21);
const juce::Colour NoiseRepellentLookAndFeel::kColorPanelBorder  = juce::Colour(0xff272c37);

NoiseRepellentLookAndFeel::NoiseRepellentLookAndFeel()
{
    setColour(juce::Slider::thumbColourId, juce::Colour(0xff303746));
    setColour(juce::Slider::trackColourId, juce::Colour(0xff0c0e12));
    setColour(juce::TextButton::textColourOffId, juce::Colour(0xff808896));
    setColour(juce::TextButton::textColourOnId, kColorNoiseProfile);
    setColour(juce::ComboBox::backgroundColourId, juce::Colour(0xff0c0e12));
    setColour(juce::ComboBox::outlineColourId, kColorPanelBorder);
}

void NoiseRepellentLookAndFeel::drawLinearSlider(juce::Graphics& g, int x, int y, int width, int height,
                                                 float sliderPos, float minSliderPos, float maxSliderPos,
                                                 const juce::Slider::SliderStyle style, juce::Slider& slider)
{
    juce::Colour fillColour = slider.findColour(juce::Slider::rotarySliderFillColourId, true);
    if (fillColour.isTransparent()) {
        fillColour = kColorFineTuning;
    }

    if (style == juce::Slider::LinearVertical)
    {
        float trackW = 4.0f;
        float trackX = x + (width - trackW) * 0.5f;

        // Background track
        g.setColour(juce::Colour(0xff090a0c));
        g.fillRoundedRectangle(trackX, (float)y, trackW, (float)height, 2.0f);

        // Solid Single-Color Fill
        g.setColour(fillColour);
        g.fillRoundedRectangle(trackX, sliderPos, trackW, (float)y + (float)height - sliderPos, 2.0f);

        // Thumb Cap
        float thumbW = 18.0f;
        float thumbH = 12.0f;
        float thumbX = x + (width - thumbW) * 0.5f;
        float thumbY = sliderPos - thumbH * 0.5f;

        g.setColour(juce::Colour(0xff303746));
        g.fillRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f);

        g.setColour(juce::Colour(0xff545f78));
        g.drawRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f, 1.0f);

        g.setColour(fillColour);
        g.fillRect(thumbX + 5.0f, thumbY + 5.0f, 8.0f, 2.0f);
    }
    else if (style == juce::Slider::LinearHorizontal)
    {
        float trackH = 4.0f;
        float trackY = y + (height - trackH) * 0.5f;

        g.setColour(juce::Colour(0xff090a0c));
        g.fillRoundedRectangle((float)x, trackY, (float)width, trackH, 2.0f);

        g.setColour(fillColour);
        g.fillRoundedRectangle((float)x, trackY, sliderPos - (float)x, trackH, 2.0f);

        float thumbW = 12.0f;
        float thumbH = 16.0f;
        float thumbX = sliderPos - thumbW * 0.5f;
        float thumbY = y + (height - thumbH) * 0.5f;

        g.setColour(juce::Colour(0xff303746));
        g.fillRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f);

        g.setColour(juce::Colour(0xff545f78));
        g.drawRoundedRectangle(thumbX, thumbY, thumbW, thumbH, 2.0f, 1.0f);

        g.setColour(fillColour);
        g.fillRect(thumbX + 5.0f, thumbY + 4.0f, 2.0f, 8.0f);
    }
}

void NoiseRepellentLookAndFeel::drawButtonBackground(juce::Graphics& g, juce::Button& button,
                                                      const juce::Colour& backgroundColour,
                                                      bool shouldDrawButtonAsHighlighted,
                                                      bool shouldDrawButtonAsDown)
{
    auto bounds = button.getLocalBounds().toFloat();
    juce::Colour base = backgroundColour;

    if (shouldDrawButtonAsDown)      base = base.darker(0.2f);
    else if (shouldDrawButtonAsHighlighted) base = base.brighter(0.1f);

    g.setColour(base);
    g.fillRoundedRectangle(bounds, 4.0f);

    g.setColour(kColorPanelBorder);
    g.drawRoundedRectangle(bounds, 4.0f, 1.0f);
}

void NoiseRepellentLookAndFeel::drawToggleButton(juce::Graphics& g, juce::ToggleButton& button,
                                                 bool shouldDrawButtonAsHighlighted,
                                                 bool shouldDrawButtonAsDown)
{
    auto bounds = button.getLocalBounds().toFloat();
    bool isOn = button.getToggleState();

    g.setColour(isOn ? kColorNoiseProfile : juce::Colour(0xff212632));
    g.fillRoundedRectangle(bounds, 4.0f);

    g.setColour(kColorPanelBorder);
    g.drawRoundedRectangle(bounds, 4.0f, 1.0f);

    g.setColour(isOn ? juce::Colour(0xff101216) : juce::Colour(0xffe6e8ed));
    g.setFont(12.0f);
    g.drawText(button.getButtonText(), bounds, juce::Justification::centred);
}
