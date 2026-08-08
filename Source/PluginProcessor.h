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

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>
#include <vector>

extern "C" {
#include "specbleach_denoiser.h"
#include "specbleach_2d_denoiser.h"
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
    void processBlockBypassed(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    juce::AudioProcessorParameter* getBypassParameter() const override;

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

    // FFT visualization constants
    static constexpr int kFftOrder = 12;                    // 2^12 = 4096 point FFT
    static constexpr size_t kFftSize = 1 << kFftOrder;      // 4096
    static constexpr size_t kFftBins = kFftSize / 2;        // 2048 unique frequency bins

    // Spectral frame shared with GUI via lock-free ring buffer
    struct SpectralFrame {
        std::array<float, kFftBins> inputMagnitudeDB{};     // dB spectrum of input
        std::array<float, kFftBins> noiseFloorDB{};         // dB spectrum of noise profile
        std::array<float, kFftBins> outputMagnitudeDB{};    // dB spectrum of output
        std::vector<float> tonalPeaksHz{};                  // Detected tonal peak frequencies in Hz
        bool hasNoiseProfile = false;
        bool isLinked = true;
    };

    bool getNextSpectralFrame(SpectralFrame& frame);

    void resetNoiseProfile();
    bool hasNoiseProfile() const;

    double getSampleRate() const { return currentSampleRate; }

private:
    void ensureEnginesInitialized(double sampleRate);
    void syncNoiseProfiles(int sourceAlgoMode);

    juce::AudioProcessorValueTreeState parameters;
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();

    // DSP Engines — one instance per channel for correct stereo processing
    SpectralBleachHandle specbleach1D_L = nullptr;
    SpectralBleachHandle specbleach1D_R = nullptr;
    SpectralBleachHandle specbleach2D_L = nullptr;
    SpectralBleachHandle specbleach2D_R = nullptr;

    juce::AudioParameterBool* bypassParameter = nullptr;
    juce::dsp::DryWetMixer<float> dryWetMixer;

    double currentSampleRate = 44100.0;
    int currentAlgoMode = 1; // Track for dynamic latency updates
    bool wasLearning = false; // Track learn mode state transition to sync profiles

    struct PendingProfile {
        int channel = 0; // 0 = Left, 1 = Right
        int mode = 1;
        uint32_t size = 0;
        uint32_t blockCount = 0;
        std::vector<float> data;
    };

    std::vector<PendingProfile> pendingProfiles;

    // Persistent dry input copy for FFT visualization (prevents RT audio thread allocation)
    std::vector<float> dryInputL;

    // Crossfading between 1D and 2D engines to prevent clicks/pops during mode changes
    juce::AudioBuffer<float> crossfadeBuffer;
    int targetAlgoMode = 1;
    float crossfadeProgress = 1.0f;
    float crossfadeStep = 0.0f;

    // FFT analysis for visualization
    juce::dsp::FFT fftAnalyzer{ kFftOrder };
    juce::dsp::WindowingFunction<float> fftWindow{ kFftSize, juce::dsp::WindowingFunction<float>::hann };
    std::array<float, kFftSize * 2> fftInputWork{};   // real+imag interleaved for input FFT
    std::array<float, kFftSize * 2> fftOutputWork{};   // real+imag interleaved for output FFT

    // Lock-free SPSC Ring Buffer for GUI visualization
    juce::AbstractFifo spectralFifo{ 16 };
    std::vector<SpectralFrame> spectralBuffer{ 16 };

    // Accumulation buffer for FFT (collects samples across processBlock calls)
    std::array<float, kFftSize> fftAccumInput{};
    std::array<float, kFftSize> fftAccumOutput{};
    size_t fftAccumCount = 0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(NoiseRepellentAudioProcessor)
};
