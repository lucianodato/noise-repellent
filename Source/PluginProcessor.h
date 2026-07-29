/*
noise-repellent -- Noise Reduction JUCE Plugin

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.
*/

#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

extern "C" {
#include "specbleach_denoiser.h"
#include "specbleach_2d_denoiser.h"
#include "DSP/signal_crossfade.h"
}

class NoiseRepellentAudioProcessor : public juce::AudioProcessor
{
public:
    NoiseRepellentAudioProcessor();
    ~NoiseRepellentAudioProcessor() override;

    void prepareToPlay(double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;

    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

    const juce::String getName() const override { return JucePlugin_Name; }

    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }
    bool isMidiEffect() const override { return false; }
    double getTailLengthSeconds() const override { return 0.0; }

    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram(int) override {}
    const juce::String getProgramName(int) override { return {}; }
    void changeProgramName(int, const juce::String&) override {}

    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;

    juce::AudioProcessorValueTreeState& getAPVTS() { return parameters; }

    // Ring buffer structure for FFT GUI visualization
    static constexpr size_t kFftSize = 512;
    struct SpectralFrame {
        std::array<float, kFftSize> inputMagnitude;
        std::array<float, kFftSize> noiseFloor;
        std::array<float, kFftSize> outputMagnitude;
        std::vector<float> tonalPeaksHz;
    };

    bool getNextSpectralFrame(SpectralFrame& frame);

    void resetNoiseProfile();

private:
    juce::AudioProcessorValueTreeState parameters;
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();

    // DSP Engines
    SpectralBleachHandle specbleach1D = nullptr;
    SpectralBleachHandle specbleach2D = nullptr;

    SignalCrossfade* softBypassL = nullptr;
    SignalCrossfade* softBypassR = nullptr;

    double currentSampleRate = 44100.0;

    // Lock-free SPSC Ring Buffer for GUI visualization
    juce::AbstractFifo spectralFifo{ 16 };
    std::vector<SpectralFrame> spectralBuffer{ 16 };

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentAudioProcessor)
};
