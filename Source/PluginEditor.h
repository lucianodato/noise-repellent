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
#include "PluginProcessor.h"
#include "GUI/LookAndFeel.h"
#include "GUI/SpectralVisualizer.h"

class NoiseRepellentAudioProcessorEditor : public juce::AudioProcessorEditor,
                                    public juce::Timer
{
public:
    explicit NoiseRepellentAudioProcessorEditor(NoiseRepellentAudioProcessor&);
    ~NoiseRepellentAudioProcessorEditor() override;

    void paint(juce::Graphics&) override;
    void resized() override;
    void timerCallback() override;

private:
    NoiseRepellentAudioProcessor& audioProcessor;
    NoiseRepellentLookAndFeel customLookAndFeel;

    // Header Controls
    juce::Label brandLabel;
    juce::ComboBox comboAlgoMode;
    juce::TextButton btnAdvancedToggle{ "Advanced Controls" };
    juce::ToggleButton btnDelta{ "Delta" };
    juce::ToggleButton btnBypass{ "Bypass" };

    // Module 1: Compact Noise Profile
    juce::GroupComponent groupProfile{ "groupProfile", "PROFILE" };
    juce::TextButton btnLearn{ "LEARN" };
    juce::TextButton btnResetProfile{ "RESET" };
    juce::ToggleButton btnAdaptiveNoise{ "ADAPTIVE" };
    juce::Label lblProfileStatus;

    // Module 2: Denoising & Spectrum
    juce::GroupComponent groupDenoising{ "groupDenoising", "DENOISING PROCESSING" };
    juce::ToggleButton btnLink{ "Linked" };
    juce::Slider sliderMasterRed{ juce::Slider::LinearVertical, juce::Slider::TextBoxBelow };
    juce::Slider sliderTonalRed{ juce::Slider::LinearVertical, juce::Slider::TextBoxBelow };
    juce::Slider sliderBroadbandSupp{ juce::Slider::LinearVertical, juce::Slider::TextBoxBelow };
    SpectralVisualizerComponent spectralVisualizer;

    // Module 3: Advanced Controls Panel (Collapsible)
    juce::GroupComponent groupAdvanced{ "groupAdvanced", "ADVANCED CONTROLS" };
    juce::ComboBox comboMethod;
    juce::Slider sliderAggressiveness{ juce::Slider::LinearHorizontal, juce::Slider::NoTextBox };
    juce::Slider sliderSmoothing{ juce::Slider::LinearHorizontal, juce::Slider::NoTextBox };
    juce::Slider sliderMasking{ juce::Slider::LinearHorizontal, juce::Slider::NoTextBox };
    juce::Slider sliderWhitening{ juce::Slider::LinearHorizontal, juce::Slider::NoTextBox };
    juce::Slider sliderSuppression{ juce::Slider::LinearHorizontal, juce::Slider::NoTextBox };

    juce::Label lblMethod{ "lblMethod", "ESTIMATION METHOD" };
    juce::Label lblAggressiveness{ "lblAggressiveness", "AGGRESSIVENESS" };
    juce::Label lblSmoothing{ "lblSmoothing", "SMOOTHING" };
    juce::Label lblMasking{ "lblMasking", "MASKING PROTECT" };
    juce::Label lblWhitening{ "lblWhitening", "WHITENING" };
    juce::Label lblSuppression{ "lblSuppression", "BROADBAND SUPP." };

    // Parameter Attachments
    using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
    using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;
    using ComboBoxAttachment = juce::AudioProcessorValueTreeState::ComboBoxAttachment;

    std::unique_ptr<ComboBoxAttachment> attachAlgoMode;
    std::unique_ptr<ButtonAttachment> attachLearn;
    std::unique_ptr<ButtonAttachment> attachAdaptive;
    std::unique_ptr<ButtonAttachment> attachLink;
    std::unique_ptr<ButtonAttachment> attachDelta;
    std::unique_ptr<ButtonAttachment> attachBypass;
    std::unique_ptr<ButtonAttachment> attachShowAdvanced;

    std::unique_ptr<ComboBoxAttachment> attachMethod;

    std::unique_ptr<SliderAttachment> attachMasterRed;
    std::unique_ptr<SliderAttachment> attachTonalRed;
    std::unique_ptr<SliderAttachment> attachBroadbandSupp;
    std::unique_ptr<SliderAttachment> attachAggressiveness;
    std::unique_ptr<SliderAttachment> attachSmoothing;
    std::unique_ptr<SliderAttachment> attachMasking;
    std::unique_ptr<SliderAttachment> attachWhitening;
    std::unique_ptr<SliderAttachment> attachSuppression;

    bool isAdvancedVisible = true;

    void updateLayout();
    void updateSliderLabels();
    void updateProfileStatus();

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentAudioProcessorEditor)
};
