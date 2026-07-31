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
    brandLabel.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeBrand, juce::Font::bold));
    brandLabel.setColour(juce::Label::textColourId, juce::Colours::white);
    addAndMakeVisible(brandLabel);

    // Algorithm Mode Label & Dropdown
    lblAlgoHeader.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
    lblAlgoHeader.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    lblAlgoHeader.setJustificationType(juce::Justification::centred);
    addAndMakeVisible(lblAlgoHeader);

    comboAlgoMode.addItemList({ "1D Spectral (Fast & Low CPU)", "2D NLM Patch (High Quality)" }, 1);
    addAndMakeVisible(comboAlgoMode);

    // Header Action Toggles
    isAdvancedVisible = audioProcessor.getAPVTS().getRawParameterValue("show_advanced")->load() > 0.5f;
    btnAdvancedToggle.setClickingTogglesState(true);
    addAndMakeVisible(btnAdvancedToggle);
    btnAdvancedToggle.setColour(juce::TextButton::buttonColourId, juce::Colour(0xff3f4757));
    btnAdvancedToggle.setColour(juce::TextButton::buttonOnColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    btnAdvancedToggle.onClick = [this]() {
        isAdvancedVisible = btnAdvancedToggle.getToggleState();
        updateLayout();
    };

    btnDelta.setClickingTogglesState(true);
    btnDelta.setColour(juce::TextButton::buttonColourId, juce::Colour(0xff3f4757));
    btnDelta.setColour(juce::TextButton::buttonOnColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    addAndMakeVisible(btnDelta);
    btnDelta.onClick = [this]() { updateLayout(); };

    btnBypass.setClickingTogglesState(true);
    btnBypass.setColour(juce::TextButton::buttonColourId, juce::Colour(0xff3f4757));
    btnBypass.setColour(juce::TextButton::buttonOnColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    addAndMakeVisible(btnBypass);
    btnBypass.onClick = [this]() { updateLayout(); };

    // Beginner Module 1: Noise Profile
    addAndMakeVisible(groupProfile);
    groupProfile.setText("NOISE PROFILE");
    groupProfile.setColour(juce::GroupComponent::outlineColourId, NoiseRepellentLookAndFeel::kColorPanelBorder);
    groupProfile.setColour(juce::GroupComponent::textColourId, NoiseRepellentLookAndFeel::kColorFineTuning);

    btnLearn.setClickingTogglesState(true);
    addAndMakeVisible(btnLearn);
    btnLearn.setColour(juce::TextButton::buttonColourId, juce::Colour(0xffc0392b));   // Prominent Studio Crimson Red CTA
    btnLearn.setColour(juce::TextButton::buttonOnColourId, juce::Colour(0xffe74c3c)); // Active Learning Red
    btnLearn.onClick = [this]() {
        bool isLearning = btnLearn.getToggleState();
        btnLearn.setButtonText(isLearning ? "Learning..." : "Learn Noise");
    };

    addAndMakeVisible(btnResetProfile);
    btnResetProfile.onClick = [this]() {
        audioProcessor.resetNoiseProfile();
    };

    addAndMakeVisible(btnAdaptiveNoise);
    btnAdaptiveNoise.onClick = [this]() {
        updateLayout();
    };

    // Profile Status Label (HUD Overlay on Chart)
    addAndMakeVisible(lblProfileStatus);
    lblProfileStatus.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
    lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    lblProfileStatus.setJustificationType(juce::Justification::left);
    lblProfileStatus.setText("STATUS: STATIONARY NOISE PROFILE", juce::dontSendNotification);

    // Beginner Module 2: Denoising & Resizable Canvas
    addAndMakeVisible(groupDenoising);
    groupDenoising.setText("");
    groupDenoising.setColour(juce::GroupComponent::outlineColourId, NoiseRepellentLookAndFeel::kColorPanelBorder);
    groupDenoising.setColour(juce::GroupComponent::textColourId, NoiseRepellentLookAndFeel::kColorFineTuning);

    lblReductionHeader.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
    lblReductionHeader.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    lblReductionHeader.setJustificationType(juce::Justification::centred);
    addAndMakeVisible(lblReductionHeader);

    addAndMakeVisible(btnLink);
    btnLink.setButtonText(btnLink.getToggleState() ? "Linked" : "Unlinked");
    btnLink.onClick = [this]() {
        bool isLinked = btnLink.getToggleState();
        btnLink.setButtonText(isLinked ? "Linked" : "Unlinked");
        if (isLinked)
        {
            // Reset reduction parameters to default (15.0 dB) when relinked
            if (auto* pMaster = audioProcessor.getAPVTS().getParameter("reduction_amount"))
                pMaster->setValueNotifyingHost(pMaster->getDefaultValue());
            if (auto* pTonal = audioProcessor.getAPVTS().getParameter("tonal_reduction"))
                pTonal->setValueNotifyingHost(pTonal->getDefaultValue());
        }
        updateLayout();
    };

    lblMasterRed.setText("BROADBAND", juce::dontSendNotification);
    lblMasterRed.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
    lblMasterRed.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    lblMasterRed.setJustificationType(juce::Justification::centred);
    addAndMakeVisible(lblMasterRed);

    lblTonalRed.setText("TONAL", juce::dontSendNotification);
    lblTonalRed.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
    lblTonalRed.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    lblTonalRed.setJustificationType(juce::Justification::centred);
    addAndMakeVisible(lblTonalRed);

    sliderMasterRed.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    sliderTonalRed.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);

    addAndMakeVisible(sliderMasterRed);
    addAndMakeVisible(sliderTonalRed);

    addAndMakeVisible(spectralVisualizer);

    // Collapsible Module 3: Advanced Controls
    addAndMakeVisible(groupAdvanced);
    groupAdvanced.setText("ADVANCED CONTROLS");
    groupAdvanced.setColour(juce::GroupComponent::outlineColourId, NoiseRepellentLookAndFeel::kColorPanelBorder);
    groupAdvanced.setColour(juce::GroupComponent::textColourId, NoiseRepellentLookAndFeel::kColorFineTuning);

    comboMethod.addItemList({ "SPP-MMSE (Unbiased)", "Brandt (Trimmed Mean)", "Martin (Min Statistics)" }, 1);
    addAndMakeVisible(comboMethod);

    lblMethod.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
    lblMethod.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    lblMethod.setJustificationType(juce::Justification::left);
    addAndMakeVisible(lblMethod);

    lblAggressiveness.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
    lblAggressiveness.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    lblAggressiveness.setJustificationType(juce::Justification::centred);
    sliderAggressiveness.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    addAndMakeVisible(sliderAggressiveness);
    addAndMakeVisible(lblAggressiveness);

    sliderSmoothing.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    sliderMasking.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    sliderWhitening.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    sliderSuppression.setColour(juce::Slider::rotarySliderFillColourId, NoiseRepellentLookAndFeel::kColorDenoising);

    addAndMakeVisible(sliderSmoothing);
    addAndMakeVisible(sliderMasking);
    addAndMakeVisible(sliderWhitening);
    addAndMakeVisible(sliderSuppression);

    for (auto* lbl : { &lblSmoothing, &lblMasking, &lblWhitening, &lblSuppression })
    {
        lbl->setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
        lbl->setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
        lbl->setJustificationType(juce::Justification::centred);
        addAndMakeVisible(lbl);
    }

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
    bool isLinked = btnLink.getToggleState();

    sliderMasterRed.setVisible(true);
    sliderTonalRed.setVisible(!isLinked);

    lblMasterRed.setVisible(!isLinked);
    lblTonalRed.setVisible(!isLinked);

    groupAdvanced.setVisible(isAdvancedVisible);
    btnAdaptiveNoise.setVisible(isAdvancedVisible);
    comboMethod.setVisible(isAdvancedVisible);
    lblMethod.setVisible(isAdvancedVisible);
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
    lblSuppression.setText("SUPPRESSION: " + juce::String(static_cast<int>(sliderSuppression.getValue())) + "%", juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::updateProfileStatus()
{
    bool isBypassed = btnBypass.getToggleState();
    bool isAdaptive = btnAdaptiveNoise.getToggleState();

    // When bypassed, disable everything except the Bypass button itself
    bool pluginActive = !isBypassed;

    // Header controls
    btnDelta.setEnabled(pluginActive);
    btnAdvancedToggle.setEnabled(pluginActive);
    comboAlgoMode.setEnabled(pluginActive);
    lblAlgoHeader.setEnabled(pluginActive);

    // Noise Profile box: disabled when bypassed OR when adaptive is active
    bool profileEnabled = pluginActive && !isAdaptive;
    groupProfile.setEnabled(pluginActive);
    btnLearn.setEnabled(profileEnabled);
    btnResetProfile.setEnabled(profileEnabled);
    sliderAggressiveness.setEnabled(profileEnabled);
    lblAggressiveness.setEnabled(profileEnabled);
    btnAdaptiveNoise.setEnabled(pluginActive);

    // Reduction controls
    btnLink.setEnabled(pluginActive);
    lblReductionHeader.setEnabled(pluginActive);
    sliderMasterRed.setEnabled(pluginActive);
    sliderTonalRed.setEnabled(pluginActive);
    lblMasterRed.setEnabled(pluginActive);
    lblTonalRed.setEnabled(pluginActive);

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

    // Status label
    bool isDelta = btnDelta.getToggleState();

    if (isBypassed) {
        lblProfileStatus.setText("STATUS: BYPASSED", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, juce::Colour(0xff808896));
    } else if (isDelta) {
        lblProfileStatus.setText("STATUS: DELTA (LISTENING TO REMOVED NOISE)", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorTonalPeaks);
    } else if (isAdaptive) {
        lblProfileStatus.setText("STATUS: ADAPTIVE NOISE ESTIMATION", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    } else {
        lblProfileStatus.setText("STATUS: STATIONARY NOISE PROFILE", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    }
}

void NoiseRepellentAudioProcessorEditor::timerCallback()
{
    updateSliderLabels();
    updateProfileStatus();
}

void NoiseRepellentAudioProcessorEditor::paint(juce::Graphics& g)
{
    g.fillAll(juce::Colour(0xff2c313c));

    if (isAdvancedVisible)
    {
        g.setColour(juce::Colour(0xff4f586c));

        float topY = (float)groupAdvanced.getY() + 28.0f;
        float bottomY = (float)groupAdvanced.getBottom() - 14.0f;

        // Separator 1: After Adaptive Noise & Estimation Method block
        if (comboMethod.getWidth() > 0)
        {
            float x1 = (float)comboMethod.getRight() + 7.0f;
            g.drawVerticalLine(juce::roundToInt(x1), topY, bottomY);
        }

        // Separator 2: After Smoothing slider block
        if (sliderSmoothing.getWidth() > 0)
        {
            float x2 = (float)sliderSmoothing.getRight() + 7.0f;
            g.drawVerticalLine(juce::roundToInt(x2), topY, bottomY);
        }

        // Separator 3: After Masking Protect slider block
        if (sliderMasking.getWidth() > 0)
        {
            float x3 = (float)sliderMasking.getRight() + 7.0f;
            g.drawVerticalLine(juce::roundToInt(x3), topY, bottomY);
        }

        // Separator 4: After Whitening slider block
        if (sliderWhitening.getWidth() > 0)
        {
            float x4 = (float)sliderWhitening.getRight() + 7.0f;
            g.drawVerticalLine(juce::roundToInt(x4), topY, bottomY);
        }
    }
}

void NoiseRepellentAudioProcessorEditor::resized()
{
    auto area = getLocalBounds().reduced(12);

    // Header (Fixed 54px with spacious headroom for profile controls)
    auto headerArea = area.removeFromTop(54);
    brandLabel.setBounds(headerArea.removeFromLeft(150).withSizeKeepingCentre(150, 26));

    btnBypass.setBounds(headerArea.removeFromRight(65).withSizeKeepingCentre(65, 24));
    headerArea.removeFromRight(6);
    btnDelta.setBounds(headerArea.removeFromRight(55).withSizeKeepingCentre(55, 24));
    headerArea.removeFromRight(12);

    // Calculate dynamic centering offset between brand logo and delta button
    constexpr int kAlgoWidth = 230;
    constexpr int kProfileWidth = 350;
    constexpr int kBoxGap = 16;
    int totalBoxesWidth = kAlgoWidth + kBoxGap + kProfileWidth;
    int availMiddleWidth = headerArea.getWidth();
    int centerOffset = std::max(0, (availMiddleWidth - totalBoxesWidth) / 2);

    headerArea.removeFromLeft(centerOffset);

    // Processing Engine Column in Header
    auto algoCol = headerArea.removeFromLeft(kAlgoWidth);
    int algoYPad = std::max(0, (algoCol.getHeight() - 44) / 2);
    algoCol.removeFromTop(algoYPad);
    lblAlgoHeader.setBounds(algoCol.removeFromTop(15));
    algoCol.removeFromTop(3);
    comboAlgoMode.setBounds(algoCol.removeFromTop(26));
    headerArea.removeFromLeft(kBoxGap);

    // Encapsulated Noise Profile Group Box in Header
    auto profileGroupArea = headerArea.removeFromLeft(kProfileWidth);
    groupProfile.setBounds(profileGroupArea);

    auto profileInner = profileGroupArea.reduced(6, 4);
    profileInner.removeFromTop(10); // Clear space for "NOISE PROFILE" group title line

    auto profileRow = profileInner;
    auto btnCol1 = profileRow.removeFromLeft(95);
    btnLearn.setBounds(btnCol1.withSizeKeepingCentre(95, 24));
    profileRow.removeFromLeft(6);

    auto btnCol2 = profileRow.removeFromLeft(55);
    btnResetProfile.setBounds(btnCol2.withSizeKeepingCentre(55, 24));
    profileRow.removeFromLeft(14);

    // Aggressiveness Slider Stack (Perfectly Vertically Centered)
    auto aggrCol = profileRow;
    constexpr int kStackHeight = 14 + 3 + 16; // 33px stack
    int aggrYPad = std::max(0, (aggrCol.getHeight() - kStackHeight) / 2);
    aggrCol.removeFromTop(aggrYPad);

    lblAggressiveness.setBounds(aggrCol.removeFromTop(14));
    aggrCol.removeFromTop(3);
    sliderAggressiveness.setBounds(aggrCol.removeFromTop(16));

    area.removeFromTop(8);

    // Bottom Collapsible Advanced Panel
    if (isAdvancedVisible)
    {
        auto advArea = area.removeFromBottom(150);
        groupAdvanced.setBounds(advArea);

        // Position Advanced Controls toggle button embedded in top-right title bar of groupAdvanced
        btnAdvancedToggle.setBounds(advArea.getRight() - 145, advArea.getY() + 2, 135, 20);

        auto advInner = advArea.reduced(12);
        advInner.removeFromTop(12);

        auto leftAdv = advInner.removeFromLeft(240);

        // Vertically center Adaptive Noise checkbox & Estimation Method controls
        constexpr int totalLeftHeight = 24 + 12 + 16 + 4 + 26; // 82px
        int vertPad = std::max(0, (leftAdv.getHeight() - totalLeftHeight) / 2);
        leftAdv.removeFromTop(vertPad);

        btnAdaptiveNoise.setBounds(leftAdv.removeFromTop(24));
        leftAdv.removeFromTop(12); // Generous padding between Adaptive checkbox and Estimation Method

        lblMethod.setBounds(leftAdv.removeFromTop(16));
        leftAdv.removeFromTop(4);
        comboMethod.setBounds(leftAdv.removeFromTop(26));

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

    // Main Denoising & Spectrum Canvas takes FULL 100% remaining width and height
    groupDenoising.setBounds(area);

    auto denoiseInner = area.reduced(10);
    denoiseInner.removeFromTop(10);

    bool isLinked = btnLink.getToggleState();
    int faderBankWidth = isLinked ? 85 : 155;

    auto faderBankArea = denoiseInner.removeFromLeft(faderBankWidth);
    lblReductionHeader.setBounds(faderBankArea.removeFromTop(16));
    faderBankArea.removeFromTop(2);
    btnLink.setBounds(faderBankArea.removeFromTop(22));
    faderBankArea.removeFromTop(6);

    if (isLinked)
    {
        lblMasterRed.setBounds(faderBankArea.removeFromTop(16));
        faderBankArea.removeFromTop(2);
        sliderMasterRed.setBounds(faderBankArea);
    }
    else
    {
        // Unlinked: show TWO sliders side by side (BROADBAND on left, TONAL on right)
        int colW = (faderBankArea.getWidth() - 6) / 2;

        auto leftCol = faderBankArea.removeFromLeft(colW);
        faderBankArea.removeFromLeft(6); // spacing between columns
        auto rightCol = faderBankArea;

        lblMasterRed.setBounds(leftCol.removeFromTop(16));
        leftCol.removeFromTop(2);
        sliderMasterRed.setBounds(leftCol);

        lblTonalRed.setBounds(rightCol.removeFromTop(16));
        rightCol.removeFromTop(2);
        sliderTonalRed.setBounds(rightCol);
    }

    denoiseInner.removeFromLeft(10);

    // Dynamic Expanding FFT Canvas takes ALL remaining full space
    spectralVisualizer.setBounds(denoiseInner);

    // Position Profile Status HUD overlay label in top-left corner of FFT display canvas (offset by 64px to clear dB Y-axis labels)
    lblProfileStatus.setBounds(spectralVisualizer.getX() + 64, spectralVisualizer.getY() + 10, 300, 20);
    lblProfileStatus.toFront(false);

    // Position Advanced Controls overlay button in top-right corner of FFT display canvas
    constexpr int advW = 148;
    constexpr int advH = 24;
    btnAdvancedToggle.setBounds(spectralVisualizer.getRight() - advW - 12, spectralVisualizer.getY() + 10, advW, advH);
    btnAdvancedToggle.toFront(false);
}
