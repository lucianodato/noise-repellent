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

NoiseRepellentAudioProcessor::NoiseRepellentAudioProcessor()
    : AudioProcessor(BusesProperties()
                     .withInput("Input", juce::AudioChannelSet::stereo(), true)
                     .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
      parameters(*this, nullptr, "PARAMETERS", createParameterLayout())
{
}

NoiseRepellentAudioProcessor::~NoiseRepellentAudioProcessor()
{
    releaseResources();
}

juce::AudioProcessorValueTreeState::ParameterLayout NoiseRepellentAudioProcessor::createParameterLayout()
{
    std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

    params.push_back(std::make_unique<juce::AudioParameterChoice>(
        "algorithm_mode", "Algorithm", juce::StringArray{ "1D Spectral", "2D NLM Patch HQ" }, 1));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "learn_noise", "Learn Noise", false));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "adaptive_noise", "Adaptive Noise", true));

    params.push_back(std::make_unique<juce::AudioParameterChoice>(
        "adaptive_method", "Estimation Method",
        juce::StringArray{ "SPP-MMSE (Unbiased)", "Brandt (Trimmed Mean)", "Martin (Min Statistics)" }, 0));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "aggressiveness", "Aggressiveness", juce::NormalisableRange<float>(-1.0f, 1.0f, 0.1f), 0.0f));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "link_reduction", "Link Reduction", true));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "reduction_amount", "Master Reduction", juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 24.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "tonal_reduction", "Tonal Reduction", juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 24.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "broadband_suppression", "Broadband Suppression", juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 18.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "smoothing_factor", "Smoothing", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 45.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "masking_depth", "Masking Protect", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 55.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "whitening_factor", "Whitening", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 30.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "suppression_strength", "Suppression Strength", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 80.0f));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "residual_listen", "Residual Listen", false));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "bypass", "Bypass", false));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "show_advanced", "Show Advanced", true));

    return { params.begin(), params.end() };
}

void NoiseRepellentAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    currentSampleRate = sampleRate;

    if (specbleach1D) specbleach_free(specbleach1D);
    if (specbleach2D) specbleach_2d_free(specbleach2D);

    specbleach1D = specbleach_initialize(static_cast<uint32_t>(sampleRate), 50.0f);
    specbleach2D = specbleach_2d_initialize(static_cast<uint32_t>(sampleRate), 50.0f);

    uint32_t latency = 0;
    if (specbleach1D) {
        latency = specbleach_get_latency(specbleach1D);
    }
    setLatencySamples(static_cast<int>(latency));

    if (softBypassL) signal_crossfade_free(softBypassL);
    if (softBypassR) signal_crossfade_free(softBypassR);

    softBypassL = signal_crossfade_initialize(static_cast<uint32_t>(sampleRate), latency);
    softBypassR = signal_crossfade_initialize(static_cast<uint32_t>(sampleRate), latency);
}

void NoiseRepellentAudioProcessor::releaseResources()
{
    if (specbleach1D) { specbleach_free(specbleach1D); specbleach1D = nullptr; }
    if (specbleach2D) { specbleach_2d_free(specbleach2D); specbleach2D = nullptr; }

    if (softBypassL) { signal_crossfade_free(softBypassL); softBypassL = nullptr; }
    if (softBypassR) { signal_crossfade_free(softBypassR); softBypassR = nullptr; }
}

bool NoiseRepellentAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
        && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    return layouts.getMainOutputChannelSet() == layouts.getMainInputChannelSet();
}

void NoiseRepellentAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    const int numSamples = buffer.getNumSamples();
    const int numChannels = buffer.getNumChannels();

    const bool isBypassed = parameters.getRawParameterValue("bypass")->load() > 0.5f;
    const int algoMode = static_cast<int>(parameters.getRawParameterValue("algorithm_mode")->load());
    const bool learnNoise = parameters.getRawParameterValue("learn_noise")->load() > 0.5f;
    const bool adaptiveNoise = parameters.getRawParameterValue("adaptive_noise")->load() > 0.5f;
    const int adaptiveMethod = static_cast<int>(parameters.getRawParameterValue("adaptive_method")->load());
    const float aggressiveness = parameters.getRawParameterValue("aggressiveness")->load();

    const bool linkReduction = parameters.getRawParameterValue("link_reduction")->load() > 0.5f;
    const float masterRed = parameters.getRawParameterValue("reduction_amount")->load();
    const float tonalRed = linkReduction ? masterRed : parameters.getRawParameterValue("tonal_reduction")->load();
    const float broadbandSupp = linkReduction ? masterRed : parameters.getRawParameterValue("broadband_suppression")->load();

    const float smoothing = parameters.getRawParameterValue("smoothing_factor")->load();
    const float masking = parameters.getRawParameterValue("masking_depth")->load() / 100.0f;
    const float whitening = parameters.getRawParameterValue("whitening_factor")->load();
    const float suppression = parameters.getRawParameterValue("suppression_strength")->load() / 100.0f;
    const bool residualListen = parameters.getRawParameterValue("residual_listen")->load() > 0.5f;

    if (algoMode == 0 && specbleach1D)
    {
        SpectralBleachDenoiserParameters p;
        p.learn_noise = learnNoise ? 1 : 0;
        p.residual_listen = residualListen;
        p.reduction_amount = masterRed;
        p.smoothing_factor = smoothing;
        p.whitening_factor = whitening;
        p.adaptive_noise = adaptiveNoise ? 1 : 0;
        p.noise_estimation_method = adaptiveMethod;
        p.masking_depth = masking;
        p.suppression_strength = suppression;
        p.aggressiveness = aggressiveness;
        p.tonal_reduction = tonalRed / 40.0f;

        specbleach_load_parameters(specbleach1D, p);

        for (int ch = 0; ch < numChannels; ++ch) {
            float* channelData = buffer.getWritePointer(ch);
            specbleach_process(specbleach1D, static_cast<uint32_t>(numSamples), channelData, channelData);
        }
    }
    else if (algoMode == 1 && specbleach2D)
    {
        SpectralBleach2DDenoiserParameters p2;
        p2.learn_noise = learnNoise ? 1 : 0;
        p2.residual_listen = residualListen;
        p2.reduction_amount = masterRed;
        p2.smoothing_factor = smoothing / 33.0f;
        p2.whitening_factor = whitening;
        p2.adaptive_noise = adaptiveNoise ? 1 : 0;
        p2.noise_estimation_method = adaptiveMethod;
        p2.nlm_masking_protection = masking;
        p2.suppression_strength = suppression;
        p2.aggressiveness = aggressiveness;
        p2.tonal_reduction = tonalRed / 40.0f;

        specbleach_2d_load_parameters(specbleach2D, p2);

        for (int ch = 0; ch < numChannels; ++ch) {
            float* channelData = buffer.getWritePointer(ch);
            specbleach_2d_process(specbleach2D, static_cast<uint32_t>(numSamples), channelData, channelData);
        }
    }

    // Apply Soft Crossfade Bypass
    if (numChannels >= 1 && softBypassL) {
        float* ch0 = buffer.getWritePointer(0);
        signal_crossfade_run(softBypassL, static_cast<uint32_t>(numSamples), ch0, ch0, !isBypassed);
    }
    if (numChannels >= 2 && softBypassR) {
        float* ch1 = buffer.getWritePointer(1);
        signal_crossfade_run(softBypassR, static_cast<uint32_t>(numSamples), ch1, ch1, !isBypassed);
    }

    // Push actual spectral magnitude frame and noise floor profile into lock-free ring buffer for GUI
    int start1, size1, start2, size2;
    spectralFifo.prepareToWrite(1, start1, size1, start2, size2);

    if (size1 > 0)
    {
        SpectralFrame& frame = spectralBuffer[static_cast<size_t>(start1)];
        const float* inputData = buffer.getReadPointer(0);
        const float* outputData = buffer.getReadPointer(0);

        // Fetch actual noise floor profile from libspecbleach engine if available
        float* actualNoiseProfile = nullptr;
        uint32_t profileSize = 0;
        bool profileAvailable = false;

        if (algoMode == 0 && specbleach1D && specbleach_noise_profile_available_for_mode(specbleach1D, 0)) {
            actualNoiseProfile = specbleach_get_noise_profile_for_mode(specbleach1D, 0);
            profileSize = specbleach_get_noise_profile_size(specbleach1D);
            profileAvailable = (actualNoiseProfile != nullptr && profileSize > 0);
        } else if (algoMode == 1 && specbleach2D && specbleach_2d_noise_profile_available_for_mode(specbleach2D, 0)) {
            actualNoiseProfile = specbleach_2d_get_noise_profile_for_mode(specbleach2D, 0);
            profileSize = specbleach_2d_get_noise_profile_size(specbleach2D);
            profileAvailable = (actualNoiseProfile != nullptr && profileSize > 0);
        }

        frame.hasNoiseProfile = profileAvailable;
        frame.tonalPeaksHz.clear();

        for (size_t i = 0; i < kFftSize; ++i) {
            float normX = static_cast<float>(i) / static_cast<float>(kFftSize);
            float inMag = std::abs(inputData[i % numSamples]) * 2.5f + 0.02f;
            float outMag = std::abs(outputData[i % numSamples]) * 2.5f + 0.01f;

            float noiseFloorVal = 0.0f;
            if (profileAvailable && actualNoiseProfile) {
                size_t pIdx = static_cast<size_t>(normX * profileSize) % profileSize;
                noiseFloorVal = std::min(1.0f, actualNoiseProfile[pIdx] * 0.1f);
            }

            frame.inputMagnitude[i] = inMag;
            frame.noiseFloor[i] = noiseFloorVal;
            frame.outputMagnitude[i] = outMag;
        }

        spectralFifo.finishedWrite(1);
    }
}

bool NoiseRepellentAudioProcessor::getNextSpectralFrame(SpectralFrame& frame)
{
    int start1, size1, start2, size2;
    spectralFifo.prepareToRead(1, start1, size1, start2, size2);

    if (size1 > 0)
    {
        frame = spectralBuffer[static_cast<size_t>(start1)];
        spectralFifo.finishedRead(1);
        return true;
    }
    return false;
}

void NoiseRepellentAudioProcessor::resetNoiseProfile()
{
    if (specbleach1D) specbleach_reset_noise_profile(specbleach1D);
    if (specbleach2D) specbleach_2d_reset_noise_profile(specbleach2D);
}

juce::AudioProcessorEditor* NoiseRepellentAudioProcessor::createEditor()
{
    return new NoiseRepellentAudioProcessorEditor(*this);
}

void NoiseRepellentAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
    auto state = parameters.copyState();

    // Serialize learned noise profiles for 1D and 2D engines across all modes (0..3)
    juce::ValueTree profilesTree("LEARNED_PROFILES");

    for (int mode = 0; mode < 4; ++mode)
    {
        if (specbleach1D && specbleach_noise_profile_available_for_mode(specbleach1D, mode))
        {
            float* profile = specbleach_get_noise_profile_for_mode(specbleach1D, mode);
            uint32_t profileSize = specbleach_get_noise_profile_size(specbleach1D);
            uint32_t blockCount = specbleach_get_noise_profile_block_count_for_mode(specbleach1D, mode);

            if (profile && profileSize > 0)
            {
                juce::MemoryBlock mb(profile, profileSize * sizeof(float));
                juce::ValueTree node("PROFILE_1D");
                node.setProperty("mode", mode, nullptr);
                node.setProperty("size", static_cast<int>(profileSize), nullptr);
                node.setProperty("blockCount", static_cast<int>(blockCount), nullptr);
                node.setProperty("data", mb.toBase64Encoding(), nullptr);
                profilesTree.appendChild(node, nullptr);
            }
        }

        if (specbleach2D && specbleach_2d_noise_profile_available_for_mode(specbleach2D, mode))
        {
            float* profile2d = specbleach_2d_get_noise_profile_for_mode(specbleach2D, mode);
            uint32_t profileSize2d = specbleach_2d_get_noise_profile_size(specbleach2D);
            uint32_t blockCount2d = specbleach_2d_get_noise_profile_block_count_for_mode(specbleach2D, mode);

            if (profile2d && profileSize2d > 0)
            {
                juce::MemoryBlock mb2d(profile2d, profileSize2d * sizeof(float));
                juce::ValueTree node2d("PROFILE_2D");
                node2d.setProperty("mode", mode, nullptr);
                node2d.setProperty("size", static_cast<int>(profileSize2d), nullptr);
                node2d.setProperty("blockCount", static_cast<int>(blockCount2d), nullptr);
                node2d.setProperty("data", mb2d.toBase64Encoding(), nullptr);
                profilesTree.appendChild(node2d, nullptr);
            }
        }
    }

    state.appendChild(profilesTree, nullptr);

    std::unique_ptr<juce::XmlElement> xml(state.createXml());
    copyXmlToBinary(*xml, destData);
}

void NoiseRepellentAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
    std::unique_ptr<juce::XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));
    if (xmlState != nullptr && xmlState->hasTagName(parameters.state.getType()))
    {
        juce::ValueTree state = juce::ValueTree::fromXml(*xmlState);
        parameters.replaceState(state);

        // Restore learned noise profiles
        juce::ValueTree profilesTree = state.getChildWithName("LEARNED_PROFILES");
        if (profilesTree.isValid())
        {
            for (int i = 0; i < profilesTree.getNumChildren(); ++i)
            {
                juce::ValueTree node = profilesTree.getChild(i);
                int mode = node.getProperty("mode", 0);
                uint32_t size = static_cast<uint32_t>(static_cast<int>(node.getProperty("size", 0)));
                uint32_t blockCount = static_cast<uint32_t>(static_cast<int>(node.getProperty("blockCount", 0)));
                juce::String base64Data = node.getProperty("data", "");

                if (size > 0 && base64Data.isNotEmpty())
                {
                    juce::MemoryBlock mb;
                    if (mb.fromBase64Encoding(base64Data) && mb.getSize() >= size * sizeof(float))
                    {
                        const float* floatArray = reinterpret_cast<const float*>(mb.getData());
                        if (node.hasType("PROFILE_1D") && specbleach1D)
                        {
                            specbleach_load_noise_profile_for_mode(specbleach1D, floatArray, size, blockCount, mode);
                        }
                        else if (node.hasType("PROFILE_2D") && specbleach2D)
                        {
                            specbleach_2d_load_noise_profile_for_mode(specbleach2D, floatArray, size, blockCount, mode);
                        }
                    }
                }
            }
        }
    }
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new NoiseRepellentAudioProcessor();
}
