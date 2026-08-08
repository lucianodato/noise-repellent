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
    setResizeLimits(840, 480, 1600, 1000);
    setSize(940, 560);

    // Brand Header
    brandLabel.setText("NOISE REPELLENT", juce::dontSendNotification);
    brandLabel.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeBrand, juce::Font::bold));
    brandLabel.setColour(juce::Label::textColourId, juce::Colours::white);
    addAndMakeVisible(brandLabel);

    // Preferences Header Button
    addAndMakeVisible(btnPreferences);
    btnPreferences.setTooltip("Preferences");
    btnPreferences.setColour(juce::TextButton::buttonColourId, juce::Colour(0xff3a4150));
    btnPreferences.onClick = [this]() {
        auto* param = audioProcessor.getAPVTS().getParameter("show_tooltips");
        if (param == nullptr)
            return;

        bool currentVal = param->getValue() > 0.5f;

        juce::PopupMenu menu;
        juce::PopupMenu::Item item("Show Tooltips");
        item.itemID = 1;
        item.isEnabled = true;
        item.isTicked = currentVal;
        item.action = [param, currentVal]() {
            param->beginChangeGesture();
            param->setValueNotifyingHost(currentVal ? 0.0f : 1.0f);
            param->endChangeGesture();
        };
        menu.addItem(item);

        menu.showMenuAsync(juce::PopupMenu::Options().withTargetComponent(btnPreferences));
    };

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
        updateProfileStatus();
    };

    addAndMakeVisible(btnResetProfile);
    btnResetProfile.onClick = [this]() {
        audioProcessor.resetNoiseProfile();
        btnLearn.setToggleState(false, juce::dontSendNotification);
        updateProfileStatus();
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

    // Footer Tooltip Bar Styling (Transparent background/border for seamless plugin background integration)
    footerTooltipLabel.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel));
    footerTooltipLabel.setColour(juce::Label::backgroundColourId, juce::Colours::transparentBlack);
    footerTooltipLabel.setColour(juce::Label::outlineColourId, juce::Colours::transparentBlack);
    footerTooltipLabel.setColour(juce::Label::textColourId, juce::Colour(0xffd0d7e2));
    footerTooltipLabel.setJustificationType(juce::Justification::centredLeft);
    footerTooltipLabel.setBorderSize(juce::BorderSize<int>(0, 4, 0, 4));
    addAndMakeVisible(footerTooltipLabel);

    // Control Tooltip Descriptions
    btnPreferences.setTooltip("Plugin preferences menu.");
    comboAlgoMode.setTooltip("Select denoising algorithm: 1D STFT (fast) or 2D NLM Patch (high quality).");
    lblAlgoHeader.setTooltip("Select denoising algorithm: 1D STFT (fast) or 2D NLM Patch (high quality).");
    btnAdvancedToggle.setTooltip("Toggle visibility of advanced DSP tuning controls.");
    btnLearn.setTooltip("Capture static noise profile from current audio input.");
    btnResetProfile.setTooltip("Clear learned noise profile.");
    btnAdaptiveNoise.setTooltip("Continuously estimate background noise without manual learning.");
    sliderAggressiveness.setTooltip("Adjust noise profile threshold offset (-1.0 to +1.0 dB).");
    lblAggressiveness.setTooltip("Adjust noise profile threshold offset (-1.0 to +1.0 dB).");
    btnDelta.setTooltip("Listen to removed noise signal (residual audio).");
    btnBypass.setTooltip("Soft bypass plugin processing with smooth wet/dry fade.");
    sliderMasterRed.setTooltip("Set total noise reduction depth in dB.");
    lblMasterRed.setTooltip("Set total noise reduction depth in dB.");
    sliderTonalRed.setTooltip("Set reduction depth specifically for tonal noise peaks.");
    lblTonalRed.setTooltip("Set reduction depth specifically for tonal noise peaks.");
    btnLink.setTooltip("Link broadband and tonal reduction controls together.");
    comboMethod.setTooltip("Select adaptive noise estimation algorithm (SPP-MMSE, Brandt, or Martin).");
    lblMethod.setTooltip("Select adaptive noise estimation algorithm (SPP-MMSE, Brandt, or Martin).");
    sliderSmoothing.setTooltip("Apply temporal smoothing across frames to reduce musical noise artifacts.");
    lblSmoothing.setTooltip("Apply temporal smoothing across frames to reduce musical noise artifacts.");
    sliderMasking.setTooltip("Adjust psychoacoustic masking threshold to protect quiet signal components.");
    lblMasking.setTooltip("Adjust psychoacoustic masking threshold to protect quiet signal components.");
    sliderWhitening.setTooltip("Adjust spectral whitening to balance residual noise coloration.");
    lblWhitening.setTooltip("Adjust spectral whitening to balance residual noise coloration.");
    sliderSuppression.setTooltip("Set maximum attenuation floor for unmasked noise frequencies.");
    lblSuppression.setTooltip("Set maximum attenuation floor for unmasked noise frequencies.");

    updateLayout();
    updateSliderLabels();
    startTimerHz(30);

    // Initial default instruction text
    auto* showTooltipsParam = audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
    if (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f)
    {
        footerTooltipLabel.setText("Loop a noise-only segment in your DAW and click Learn Noise to capture a profile, then adjust reduction.", juce::dontSendNotification);
    }

    // Register parameter listener for show_tooltips
    audioProcessor.getAPVTS().addParameterListener("show_tooltips", this);

    // Register mouse listener AFTER all child components are added
    addMouseListener(this, true);
}

void NoiseRepellentAudioProcessorEditor::mouseEnter(const juce::MouseEvent& event)
{
    auto* showTooltipsParam = audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
    if (showTooltipsParam != nullptr && showTooltipsParam->load() < 0.5f)
    {
        footerTooltipLabel.setText({}, juce::dontSendNotification);
        return;
    }

    auto* comp = event.originalComponent;
    while (comp != nullptr)
    {
        if (auto* client = dynamic_cast<juce::TooltipClient*>(comp))
        {
            auto tip = client->getTooltip();
            if (tip.isNotEmpty())
            {
                footerTooltipLabel.setText(tip, juce::dontSendNotification);
                return;
            }
        }
        comp = comp->getParentComponent();
    }

    // Default instruction when hovering over background / non-tooltip components
    footerTooltipLabel.setText("Loop a noise-only segment in your DAW and click Learn Noise to capture a profile, then adjust reduction.", juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::mouseMove(const juce::MouseEvent& event)
{
    mouseEnter(event);
}

void NoiseRepellentAudioProcessorEditor::mouseExit(const juce::MouseEvent& /*event*/)
{
    auto* showTooltipsParam = audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
    if (showTooltipsParam != nullptr && showTooltipsParam->load() > 0.5f)
    {
        footerTooltipLabel.setText("Loop a noise-only segment in your DAW and click Learn Noise to capture a profile, then adjust reduction.", juce::dontSendNotification);
    }
    else
    {
        footerTooltipLabel.setText({}, juce::dontSendNotification);
    }
}

NoiseRepellentAudioProcessorEditor::~NoiseRepellentAudioProcessorEditor()
{
    audioProcessor.getAPVTS().removeParameterListener("show_tooltips", this);
    cancelPendingUpdate();
    stopTimer();
    setLookAndFeel(nullptr);
}

void NoiseRepellentAudioProcessorEditor::parameterChanged(const juce::String& parameterID, float /*newValue*/)
{
    if (parameterID == "show_tooltips")
    {
        tooltipStateDirty.store(true, std::memory_order_relaxed);
    }
}

void NoiseRepellentAudioProcessorEditor::handleAsyncUpdate()
{
    auto* showTooltipsParam = audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
    if (showTooltipsParam != nullptr && showTooltipsParam->load() > 0.5f)
    {
        footerTooltipLabel.setText("Loop a noise-only segment in your DAW and click Learn Noise to capture a profile, then adjust reduction.", juce::dontSendNotification);
    }
    else
    {
        footerTooltipLabel.setText({}, juce::dontSendNotification);
    }
    updateLayout();
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
    bool isCompact = getWidth() < 900;
    lblAggressiveness.setText(isCompact ? "AGGR: " + juce::String(sliderAggressiveness.getValue(), 1) : "AGGRESSIVENESS: " + juce::String(sliderAggressiveness.getValue(), 1), juce::dontSendNotification);
    lblSmoothing.setText("SMOOTHING: " + juce::String(static_cast<int>(sliderSmoothing.getValue())) + "%", juce::dontSendNotification);
    lblMasking.setText(isCompact ? "MASKING: " + juce::String(static_cast<int>(sliderMasking.getValue())) + "%" : "MASKING PROTECT: " + juce::String(static_cast<int>(sliderMasking.getValue())) + "%", juce::dontSendNotification);
    lblWhitening.setText("WHITENING: " + juce::String(static_cast<int>(sliderWhitening.getValue())) + "%", juce::dontSendNotification);
    lblSuppression.setText("SUPPRESSION: " + juce::String(static_cast<int>(sliderSuppression.getValue())) + "%", juce::dontSendNotification);
}

void NoiseRepellentAudioProcessorEditor::updateProfileStatus()
{
    bool isBypassed = btnBypass.getToggleState();
    bool isAdaptive = btnAdaptiveNoise.getToggleState();
    bool isLearning = btnLearn.getToggleState();
    bool hasProfile = audioProcessor.hasNoiseProfile();

    if (isLearning)
    {
        btnLearn.setButtonText("Learning...");
        btnLearn.setColour(juce::TextButton::buttonColourId, juce::Colour(0xffe74c3c));   // Active Learning Red
        btnLearn.setColour(juce::TextButton::buttonOnColourId, juce::Colour(0xffe74c3c));
    }
    else if (hasProfile)
    {
        // Profile is active: match golden amber yellow profile curve & set text to "Re-Learn"
        btnLearn.setButtonText("Re-Learn");
        btnLearn.setColour(juce::TextButton::buttonColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
        btnLearn.setColour(juce::TextButton::buttonOnColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    }
    else
    {
        // No profile: prominent red CTA "Learn Noise"
        btnLearn.setButtonText("Learn Noise");
        btnLearn.setColour(juce::TextButton::buttonColourId, juce::Colour(0xffc0392b));   // Prominent Crimson Red CTA
        btnLearn.setColour(juce::TextButton::buttonOnColourId, juce::Colour(0xffe74c3c)); // Active Learning Red
    }
    btnLearn.repaint();

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
        lblProfileStatus.setText("STATUS: DELTA", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorTonalPeaks);
    } else if (isAdaptive) {
        lblProfileStatus.setText("STATUS: ADAPTIVE", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorDenoising);
    } else if (hasProfile) {
        lblProfileStatus.setText("STATUS: PROFILE ACTIVE", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, NoiseRepellentLookAndFeel::kColorNoiseProfile);
    } else {
        lblProfileStatus.setText("STATUS: NO PROFILE", juce::dontSendNotification);
        lblProfileStatus.setColour(juce::Label::textColourId, juce::Colour(0xff808896));
    }
}

void NoiseRepellentAudioProcessorEditor::timerCallback()
{
    if (tooltipStateDirty.exchange(false, std::memory_order_relaxed))
    {
        triggerAsyncUpdate();
    }

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
        if (comboMethod.getWidth() > 0 && sliderSmoothing.getWidth() > 0)
        {
            float x1 = (float)(comboMethod.getRight() + sliderSmoothing.getX()) / 2.0f;
            g.drawVerticalLine(juce::roundToInt(x1), topY, bottomY);
        }

        // Separator 2: After Smoothing slider block
        if (sliderSmoothing.getWidth() > 0 && sliderMasking.getWidth() > 0)
        {
            float x2 = (float)(sliderSmoothing.getRight() + sliderMasking.getX()) / 2.0f;
            g.drawVerticalLine(juce::roundToInt(x2), topY, bottomY);
        }

        // Separator 3: After Masking Protect slider block
        if (sliderMasking.getWidth() > 0 && sliderWhitening.getWidth() > 0)
        {
            float x3 = (float)(sliderMasking.getRight() + sliderWhitening.getX()) / 2.0f;
            g.drawVerticalLine(juce::roundToInt(x3), topY, bottomY);
        }

        // Separator 4: After Whitening slider block
        if (sliderWhitening.getWidth() > 0 && sliderSuppression.getWidth() > 0)
        {
            float x4 = (float)(sliderWhitening.getRight() + sliderSuppression.getX()) / 2.0f;
            g.drawVerticalLine(juce::roundToInt(x4), topY, bottomY);
        }
    }
}

void NoiseRepellentAudioProcessorEditor::resized()
{
    auto area = getLocalBounds().reduced(12);

    // Header (Fixed 54px with spacious headroom for profile controls)
    auto headerArea = area.removeFromTop(54);

    // Left Title & Options Block
    brandLabel.setBounds(headerArea.removeFromLeft(135).withSizeKeepingCentre(135, 26));
    headerArea.removeFromLeft(4);
    btnPreferences.setBounds(headerArea.removeFromLeft(22).withSizeKeepingCentre(22, 22));

    // Right Action Buttons Block
    btnBypass.setBounds(headerArea.removeFromRight(60).withSizeKeepingCentre(60, 24));
    headerArea.removeFromRight(6);
    btnDelta.setBounds(headerArea.removeFromRight(50).withSizeKeepingCentre(50, 24));

    // Forced consistent horizontal gap between all header modules
    constexpr int kHeaderGap = 16;

    // 16px gap between options button and algorithm dropdown
    headerArea.removeFromLeft(kHeaderGap);

    // Calculate middle module widths so gaps on left, middle, and right are all identically kHeaderGap
    int availMiddleWidth = headerArea.getWidth(); // Width remaining up to delta button
    int totalModuleWidth = availMiddleWidth - kHeaderGap; // Available width for algo column + kHeaderGap + profile box

    int kAlgoWidth = std::clamp((totalModuleWidth - kHeaderGap) * 38 / 100, 170, 220);
    int kProfileWidth = totalModuleWidth - kHeaderGap - kAlgoWidth;

    // Processing Engine Column in Header
    auto algoCol = headerArea.removeFromLeft(kAlgoWidth);
    int algoYPad = std::max(0, (algoCol.getHeight() - 44) / 2);
    algoCol.removeFromTop(algoYPad);
    lblAlgoHeader.setBounds(algoCol.removeFromTop(15));
    algoCol.removeFromTop(3);
    comboAlgoMode.setBounds(algoCol.removeFromTop(26));

    // 16px gap between algorithm dropdown and noise profile box
    headerArea.removeFromLeft(kHeaderGap);

    // Encapsulated Noise Profile Group Box in Header (ends with 16px gap before delta button)
    auto profileGroupArea = headerArea.removeFromLeft(kProfileWidth);
    groupProfile.setBounds(profileGroupArea);

    auto profileInner = profileGroupArea.reduced(6, 4);
    profileInner.removeFromTop(10); // Clear space for "NOISE PROFILE" group title line

    auto profileRow = profileInner;
    int learnW = std::clamp(profileGroupArea.getWidth() * 28 / 100, 75, 95);
    auto btnCol1 = profileRow.removeFromLeft(learnW);
    btnLearn.setBounds(btnCol1.withSizeKeepingCentre(learnW, 24));
    profileRow.removeFromLeft(4);

    int resetW = std::clamp(profileGroupArea.getWidth() * 18 / 100, 48, 55);
    auto btnCol2 = profileRow.removeFromLeft(resetW);
    btnResetProfile.setBounds(btnCol2.withSizeKeepingCentre(resetW, 24));
    profileRow.removeFromLeft(6);

    // Aggressiveness Slider Stack (Perfectly Vertically Centered)
    auto aggrCol = profileRow;
    constexpr int kStackHeight = 14 + 3 + 16; // 33px stack
    int aggrYPad = std::max(0, (aggrCol.getHeight() - kStackHeight) / 2);
    aggrCol.removeFromTop(aggrYPad);

    lblAggressiveness.setBounds(aggrCol.removeFromTop(14));
    aggrCol.removeFromTop(3);
    sliderAggressiveness.setBounds(aggrCol.removeFromTop(16));

    area.removeFromTop(8);

    // Carve footer bar off the very bottom of outer area only when tooltips preference is enabled
    auto* showTooltipsParam = audioProcessor.getAPVTS().getRawParameterValue("show_tooltips");
    bool showTooltips = (showTooltipsParam == nullptr || showTooltipsParam->load() > 0.5f);

    footerTooltipLabel.setVisible(showTooltips);

    if (showTooltips)
    {
        constexpr int kFooterHeight = 24;
        auto footerArea = area.removeFromBottom(kFooterHeight);
        area.removeFromBottom(8);
        footerTooltipLabel.setBounds(footerArea);
    }

    // Bottom Collapsible Advanced Panel
    if (isAdvancedVisible)
    {
        auto advArea = area.removeFromBottom(150);
        groupAdvanced.setBounds(advArea);

        // Position Advanced Controls toggle button embedded in top-right title bar of groupAdvanced
        btnAdvancedToggle.setBounds(advArea.getRight() - 145, advArea.getY() + 2, 135, 20);

        auto advInner = advArea.reduced(12);
        advInner.removeFromTop(12);

        // Dynamic responsive column width allocation
        constexpr int kNumGaps = 4;
        constexpr int kGapW = 12;
        int totalAvailWidth = advInner.getWidth() - (kNumGaps * kGapW);

        // Allocate ~26% of available width to Left Block (Adaptive & Method), clamped between 160px and 220px
        int leftAdvW = std::clamp(totalAvailWidth * 26 / 100, 160, 220);

        // Distribute remaining width equally among the 4 slider columns
        int sliderW = (totalAvailWidth - leftAdvW) / 4;

        auto leftAdv = advInner.removeFromLeft(leftAdvW);

        // Vertically center Adaptive Noise checkbox & Estimation Method controls
        constexpr int totalLeftHeight = 24 + 12 + 16 + 4 + 26; // 82px
        int vertPad = std::max(0, (leftAdv.getHeight() - totalLeftHeight) / 2);
        leftAdv.removeFromTop(vertPad);

        btnAdaptiveNoise.setBounds(leftAdv.removeFromTop(24));
        leftAdv.removeFromTop(12);

        lblMethod.setBounds(leftAdv.removeFromTop(16));
        leftAdv.removeFromTop(4);
        comboMethod.setBounds(leftAdv.removeFromTop(26));

        advInner.removeFromLeft(kGapW);

        auto s1 = advInner.removeFromLeft(sliderW);
        lblSmoothing.setBounds(s1.removeFromTop(18));
        s1.removeFromTop(4);
        sliderSmoothing.setBounds(s1);

        advInner.removeFromLeft(kGapW);
        auto s2 = advInner.removeFromLeft(sliderW);
        lblMasking.setBounds(s2.removeFromTop(18));
        s2.removeFromTop(4);
        sliderMasking.setBounds(s2);

        advInner.removeFromLeft(kGapW);
        auto s3 = advInner.removeFromLeft(sliderW);
        lblWhitening.setBounds(s3.removeFromTop(18));
        s3.removeFromTop(4);
        sliderWhitening.setBounds(s3);

        advInner.removeFromLeft(kGapW);
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

    // Position Profile Status HUD overlay label in top-left corner of FFT display canvas (offset by 50px to clear dB Y-axis labels)
    lblProfileStatus.setBounds(spectralVisualizer.getX() + 50, spectralVisualizer.getY() + 10, 140, 20);
    lblProfileStatus.toFront(false);

    // Position Advanced Controls overlay button in top-right corner of FFT display canvas
    constexpr int advW = 148;
    constexpr int advH = 24;
    btnAdvancedToggle.setBounds(spectralVisualizer.getRight() - advW - 12, spectralVisualizer.getY() + 10, advW, advH);
    btnAdvancedToggle.toFront(false);

    btnPreferences.toFront(false);
    footerTooltipLabel.toFront(false);
}
