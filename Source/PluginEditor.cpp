/*
noise-repellent -- Noise Reduction JUCE Plugin

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

NoiseRepellentAudioProcessorEditor::NoiseRepellentAudioProcessorEditor(NoiseRepellentAudioProcessor& p)
    : AudioProcessorEditor(&p), audioProcessor(p), spectralVisualizer(p)
{
    setLookAndFeel(&customLookAndFeel);

    // 1. Native Resizable Window Settings
    setResizable(true, true);
    setResizeLimits(800, 480, 1400, 900);
    setSize(980, 580);

    // Brand Header
    brandLabel.setText("NOISE REPELLENT", juce::dontSendNotification);
    brandLabel.setFont(juce::FontOptions(18.0f, juce::Font::bold));
    brandLabel.setColour(juce::Label::textColourId, juce::Colours::white);
    addAndMakeVisible(brandLabel);

    // Algorithm Mode Dropdown (with Quality indicators)
    comboAlgoMode.addItemList({ "⚡ 1D Spectral (Fast)", "✨ 2D NLM Patch (High Quality)" }, 1);
    addAndMakeVisible(comboAlgoMode);

    // Header Action Toggles
    isAdvancedVisible = audioProcessor.getAPVTS().getRawParameterValue("show_advanced")->load() > 0.5f;
    btnAdvancedToggle.setClickingTogglesState(true);
    addAndMakeVisible(btnAdvancedToggle);
    btnAdvancedToggle.onClick = [this]() {
        isAdvancedVisible = btnAdvancedToggle.getToggleState();
        updateLayout();
    };

    addAndMakeVisible(btnDelta);
    addAndMakeVisible(btnBypass);

    // Beginner Module 1: Noise Profile
    addAndMakeVisible(groupProfile);
    groupProfile.setText("NOISE PROFILE");
    groupProfile.setColour(juce::GroupComponent::outlineColourId, NoiseRepellentLookAndFeel::kColorPanelBorder);
    groupProfile.setColour(juce::GroupComponent::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);

    btnLearn.setClickingTogglesState(true);
    addAndMakeVisible(btnLearn);
    btnLearn.setColour(juce::TextButton::buttonColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btnLearn.setColour(juce::TextButton::buttonOnColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile.brighter(0.25f));
    btnLearn.onClick = [this]() {
        bool isLearning = btnLearn.getToggleState();
        btnLearn.setButtonText(isLearning ? "LEARNING..." : "LEARN NOISE");
    };

    addAndMakeVisible(btnResetProfile);
    btnResetProfile.onClick = [this]() {
        audioProcessor.resetNoiseProfile();
    };

    addAndMakeVisible(btnAdaptiveNoise);

    // Profile Status Info Card (Utilizing vertical space)
    addAndMakeVisible(lblProfileStatus);
    lblProfileStatus.setFont(juce::FontOptions(11.0f, juce::Font::bold));
    lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    lblProfileStatus.setText("STATUS: NO PROFILE", juce::dontSendNotification);

    addAndMakeVisible(lblProfileInfo);
    lblProfileInfo.setFont(juce::FontOptions(10.0f));
    lblProfileInfo.setColour(juce::Label::textColourId, juce::Colour(0xff94a3b8));
    lblProfileInfo.setText("Click 'LEARN NOISE' to capture stationary noise", juce::dontSendNotification);

    // Beginner Module 2: Denoising & Resizable Canvas
    addAndMakeVisible(groupDenoising);
    groupDenoising.setText("DENOISING PROCESSING");
    groupDenoising.setColour(juce::GroupComponent::outlineColourId, NoiseRepellentLookAndFeel::kColorPanelBorder);
    groupDenoising.setColour(juce::GroupComponent::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);

    addAndMakeVisible(btnLink);
    btnLink.onClick = [this]() {
        bool isLinked = btnLink.getToggleState();
        sliderTonalRed.setVisible(!isLinked);
        sliderBroadbandSupp.setVisible(!isLinked);
        sliderMasterRed.setVisible(isLinked);
        resized();
    };

    sliderMasterRed.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    sliderTonalRed.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    sliderBroadbandSupp.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);

    addAndMakeVisible(sliderMasterRed);
    addAndMakeVisible(sliderTonalRed);
    addAndMakeVisible(sliderBroadbandSupp);

    sliderTonalRed.setVisible(false);
    sliderBroadbandSupp.setVisible(false);

    addAndMakeVisible(spectralVisualizer);

    // Collapsible Module 3: Advanced Controls
    addAndMakeVisible(groupAdvanced);
    groupAdvanced.setText("ADVANCED CONTROLS");
    groupAdvanced.setColour(juce::GroupComponent::outlineColourId, NoiseRepellentLookAndFeel::kColorPanelBorder);
    groupAdvanced.setColour(juce::GroupComponent::textColourId, NoiseRepellentLookAndFeel::kColorFineTuning);

    comboMethod.addItemList({ "SPP-MMSE (Unbiased)", "Brandt (Trimmed Mean)", "Martin (Min Statistics)" }, 1);
    addAndMakeVisible(comboMethod);
    addAndMakeVisible(lblMethod);

    sliderAggressiveness.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    addAndMakeVisible(sliderAggressiveness);
    addAndMakeVisible(lblAggressiveness);

    sliderSmoothing.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorFineTuning);
    sliderMasking.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorFineTuning);
    sliderWhitening.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorFineTuning);
    sliderSuppression.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorFineTuning);

    addAndMakeVisible(sliderSmoothing);
    addAndMakeVisible(sliderMasking);
    addAndMakeVisible(sliderWhitening);
    addAndMakeVisible(sliderSuppression);

    addAndMakeVisible(lblSmoothing);
    addAndMakeVisible(lblMasking);
    addAndMakeVisible(lblWhitening);
    addAndMakeVisible(lblSuppression);

    // Parameter APVTS Attachments
    auto& apvts = audioProcessor.getAPVTS();

    attachAlgoMode = std::make_unique<ComboBoxAttachment>(apvts, "algorithm_mode", comboAlgoMode);
    attachLearn = std::make_unique<ButtonAttachment>(apvts, "learn_noise", btnLearn);
    attachAdaptive = std::make_unique<ButtonAttachment>(apvts, "adaptive_noise", btnAdaptiveNoise);
    attachLink = std::make_unique<ButtonAttachment>(apvts, "link_reduction", btnLink);
    attachDelta = std::make_unique<ButtonAttachment>(apvts, "residual_listen", btnDelta);
    attachBypass = std::make_unique<ButtonAttachment>(apvts, "bypass", btnBypass);
    attachShowAdvanced = std::make_unique<ButtonAttachment>(apvts, "show_advanced", btnAdvancedToggle);

    attachMethod = std::make_unique<ComboBoxAttachment>(apvts, "adaptive_method", comboMethod);

    attachMasterRed = std::make_unique<SliderAttachment>(apvts, "reduction_amount", sliderMasterRed);
    attachTonalRed = std::make_unique<SliderAttachment>(apvts, "tonal_reduction", sliderTonalRed);
    attachBroadbandSupp = std::make_unique<SliderAttachment>(apvts, "broadband_suppression", sliderBroadbandSupp);
    attachAggressiveness = std::make_unique<SliderAttachment>(apvts, "aggressiveness", sliderAggressiveness);
    attachSmoothing = std::make_unique<SliderAttachment>(apvts, "smoothing_factor", sliderSmoothing);
    attachMasking = std::make_unique<SliderAttachment>(apvts, "masking_depth", sliderMasking);
    attachWhitening = std::make_unique<SliderAttachment>(apvts, "whitening_factor", sliderWhitening);
    attachSuppression = std::make_unique<SliderAttachment>(apvts, "suppression_strength", sliderSuppression);

    updateLayout();
    updateSliderLabels();
    startTimerHz(30);
}

NoiseRepellentAudioProcessorEditor::~NoiseRepellentAudioProcessorEditor()
{
    stopTimer();
    setLookAndFeel(nullptr);
}

void NoiseRepellentAudioProcessorEditor::updateLayout()
{
    groupAdvanced.setVisible(isAdvancedVisible);
    comboMethod.setVisible(isAdvancedVisible);
    lblMethod.setVisible(isAdvancedVisible);
    sliderAggressiveness.setVisible(isAdvancedVisible);
    lblAggressiveness.setVisible(isAdvancedVisible);
    sliderSmoothing.setVisible(isAdvancedVisible);
    lblSmoothing.setVisible(isAdvancedVisible);
    sliderMasking.setVisible(isAdvancedVisible);
    lblMasking.setVisible(isAdvancedVisible);
    sliderWhitening.setVisible(isAdvancedVisible);
    lblWhitening.setVisible(isAdvancedVisible);
    sliderSuppression.setVisible(isAdvancedVisible);
    lblSuppression.setVisible(isAdvancedVisible);

    resized();
}

void NoiseRepellentAudioProcessorEditor::updateSliderLabels()
{
    lblAggressiveness.setText("AGGRESSIVENESS: " + juce::String(sliderAggressiveness.getValue(), 1), juce::dontSendNotification);
    lblSmoothing.setText("SMOOTHING: " + juce::String(static_cast<int>(sliderSmoothing.getValue())) + "%", juce::dontSendNotification);
    lblMasking.setText("MASKING PROTECT: " + juce::String(static_cast<int>(sliderMasking.getValue())) + "%", juce::dontSendNotification);
    lblWhitening.setText("WHITENING: " + juce::String(static_cast<int>(sliderWhitening.getValue())) + "%", juce::dontSendNotification);
    lblSuppression.setText("BROADBAND SUPP.: " + juce::String(static_cast<int>(sliderSuppression.getValue())) + "%", juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::updateProfileStatus()
{
    bool isAdaptive = btnAdaptiveNoise.getToggleState();

    if (isAdaptive) {
        lblProfileStatus.setText("STATUS: ADAPTIVE ESTIMATION", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
        lblProfileInfo.setText("Noise floor estimated continuously from audio input", juce::dontSendNotification);
    } else {
        lblProfileStatus.setText("STATUS: STATIONARY PROFILE", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
        lblProfileInfo.setText("Using stationary noise profile snapshot", juce::dontSendNotification);
    }
}

void NoiseRepellentAudioProcessorEditor::timerCallback()
{
    updateSliderLabels();
    updateProfileStatus();
}

void NoiseRepellentAudioProcessorEditor::paint(juce::Graphics& g)
{
    g.fillAll(juce::Colour(0xff12141a));
}

void NoiseRepellentAudioProcessorEditor::resized()
{
    auto area = getLocalBounds().reduced(12);

    // Header (Fixed 36px)
    auto headerArea = area.removeFromTop(36);
    brandLabel.setBounds(headerArea.removeFromLeft(170));

    btnBypass.setBounds(headerArea.removeFromRight(65));
    headerArea.removeFromRight(6);
    btnDelta.setBounds(headerArea.removeFromRight(55));
    headerArea.removeFromRight(10);
    btnAdvancedToggle.setBounds(headerArea.removeFromRight(145));

    comboAlgoMode.setBounds(headerArea.removeFromLeft(220));

    area.removeFromTop(8);

    // Bottom Collapsible Advanced Panel
    if (isAdvancedVisible)
    {
        auto advArea = area.removeFromBottom(150);
        groupAdvanced.setBounds(advArea);

        auto advInner = advArea.reduced(12);
        advInner.removeFromTop(12);

        auto leftAdv = advInner.removeFromLeft(240);
        lblMethod.setBounds(leftAdv.removeFromTop(18));
        leftAdv.removeFromTop(4);
        comboMethod.setBounds(leftAdv.removeFromTop(26));
        leftAdv.removeFromTop(8);
        lblAggressiveness.setBounds(leftAdv.removeFromTop(18));
        leftAdv.removeFromTop(4);
        sliderAggressiveness.setBounds(leftAdv.removeFromTop(26));

        advInner.removeFromLeft(14);

        int sliderW = (advInner.getWidth() - 42) / 4;

        auto s1 = advInner.removeFromLeft(sliderW);
        lblSmoothing.setBounds(s1.removeFromTop(18));
        s1.removeFromTop(4);
        sliderSmoothing.setBounds(s1);

        advInner.removeFromLeft(14);
        auto s2 = advInner.removeFromLeft(sliderW);
        lblMasking.setBounds(s2.removeFromTop(18));
        s2.removeFromTop(4);
        sliderMasking.setBounds(s2);

        advInner.removeFromLeft(14);
        auto s3 = advInner.removeFromLeft(sliderW);
        lblWhitening.setBounds(s3.removeFromTop(18));
        s3.removeFromTop(4);
        sliderWhitening.setBounds(s3);

        advInner.removeFromLeft(14);
        lblSuppression.setBounds(advInner.removeFromTop(18));
        advInner.removeFromTop(4);
        sliderSuppression.setBounds(advInner);

        area.removeFromBottom(8);
    }

    // Beginner Main Grid (Profile + Denoising & Resizable Canvas)
    auto profileArea = area.removeFromLeft(220);
    groupProfile.setBounds(profileArea);

    auto profileInner = profileArea.reduced(12);
    profileInner.removeFromTop(16);

    btnLearn.setBounds(profileInner.removeFromTop(36));
    profileInner.removeFromTop(8);
    btnResetProfile.setBounds(profileInner.removeFromTop(30));
    profileInner.removeFromTop(8);
    btnAdaptiveNoise.setBounds(profileInner.removeFromTop(28));
    profileInner.removeFromTop(14);

    lblProfileStatus.setBounds(profileInner.removeFromTop(18));
    profileInner.removeFromTop(2);
    lblProfileInfo.setBounds(profileInner.removeFromTop(40));

    area.removeFromLeft(12);
    groupDenoising.setBounds(area);

    auto denoiseInner = area.reduced(10);
    denoiseInner.removeFromTop(16);

    auto faderBankArea = denoiseInner.removeFromLeft(85);
    btnLink.setBounds(faderBankArea.removeFromTop(24));
    faderBankArea.removeFromTop(6);

    if (sliderMasterRed.isVisible())
    {
        sliderMasterRed.setBounds(faderBankArea);
    }
    else
    {
        int faderW = (faderBankArea.getWidth() - 4) / 2;
        sliderTonalRed.setBounds(faderBankArea.removeFromLeft(faderW));
        faderBankArea.removeFromLeft(4);
        sliderBroadbandSupp.setBounds(faderBankArea);
    }

    denoiseInner.removeFromLeft(10);

    // Dynamic Expanding FFT Canvas takes ALL remaining space
    spectralVisualizer.setBounds(denoiseInner);
}
