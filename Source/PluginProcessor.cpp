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
#include <cmath>
#include <cstring>

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
        "adaptive_noise", "Adaptive Noise", false));

    params.push_back(std::make_unique<juce::AudioParameterChoice>(
        "adaptive_method", "Estimation Method",
        juce::StringArray{ "SPP-MMSE (Unbiased)", "Brandt (Trimmed Mean)", "Martin (Min Statistics)" }, 2));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "aggressiveness", "Aggressiveness", juce::NormalisableRange<float>(-1.0f, 1.0f, 0.1f), 0.0f));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "link_reduction", "Link Reduction", true));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "reduction_amount", "Master Reduction", juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 15.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "tonal_reduction", "Tonal Reduction", juce::NormalisableRange<float>(0.0f, 40.0f, 0.1f), 15.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "smoothing_factor", "Smoothing", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 0.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "nlm_smoothing", "NLM Smoothing", juce::NormalisableRange<float>(0.1f, 10.0f, 0.1f), 5.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "masking_depth", "Masking Protect", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 100.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "whitening_factor", "Whitening", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 0.0f));

    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        "suppression_strength", "Suppression", juce::NormalisableRange<float>(0.0f, 100.0f, 1.0f), 50.0f));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "residual_listen", "Residual Listen", false));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "bypass", "Bypass", false));

    params.push_back(std::make_unique<juce::AudioParameterBool>(
        "show_advanced", "Show Advanced Controls", false));

    return { params.begin(), params.end() };
}

void NoiseRepellentAudioProcessor::ensureEnginesInitialized(double sampleRate)
{
    bool needNewEngines = (specbleach1D_L == nullptr || std::abs(currentSampleRate - sampleRate) > 0.001);

    if (!needNewEngines)
        return;

    struct SavedProfileData {
        int mode;
        uint32_t size;
        uint32_t blockCount;
        std::vector<float> data;
    };

    std::vector<SavedProfileData> profiles1D;
    std::vector<SavedProfileData> profiles2D;

    // Backup 1D profiles if handles exist
    if (specbleach1D_L)
    {
        for (int mode = 1; mode <= 4; ++mode)
        {
            if (specbleach_noise_profile_available_for_mode(specbleach1D_L, mode))
            {
                float* p = specbleach_get_noise_profile_for_mode(specbleach1D_L, mode);
                uint32_t sz = specbleach_get_noise_profile_size(specbleach1D_L);
                uint32_t bc = specbleach_get_noise_profile_block_count_for_mode(specbleach1D_L, mode);
                if (p && sz > 0)
                {
                    profiles1D.push_back({ mode, sz, bc, std::vector<float>(p, p + sz) });
                }
            }
        }
    }

    // Backup 2D profiles if handles exist
    if (specbleach2D_L)
    {
        for (int mode = 1; mode <= 4; ++mode)
        {
            if (specbleach_2d_noise_profile_available_for_mode(specbleach2D_L, mode))
            {
                float* p = specbleach_2d_get_noise_profile_for_mode(specbleach2D_L, mode);
                uint32_t sz = specbleach_2d_get_noise_profile_size(specbleach2D_L);
                uint32_t bc = specbleach_2d_get_noise_profile_block_count_for_mode(specbleach2D_L, mode);
                if (p && sz > 0)
                {
                    profiles2D.push_back({ mode, sz, bc, std::vector<float>(p, p + sz) });
                }
            }
        }
    }

    // Free existing handles
    if (specbleach1D_L) specbleach_free(specbleach1D_L);
    if (specbleach1D_R) specbleach_free(specbleach1D_R);
    if (specbleach2D_L) specbleach_2d_free(specbleach2D_L);
    if (specbleach2D_R) specbleach_2d_free(specbleach2D_R);

    // Initialize per-channel instances
    specbleach1D_L = specbleach_initialize(static_cast<uint32_t>(sampleRate), 50.0f);
    specbleach1D_R = specbleach_initialize(static_cast<uint32_t>(sampleRate), 50.0f);
    specbleach2D_L = specbleach_2d_initialize(static_cast<uint32_t>(sampleRate), 50.0f);
    specbleach2D_R = specbleach_2d_initialize(static_cast<uint32_t>(sampleRate), 50.0f);

    currentSampleRate = sampleRate;

    // Restore saved profiles into new handles
    for (const auto& item : profiles1D)
    {
        if (specbleach1D_L)
            specbleach_load_noise_profile_for_mode(specbleach1D_L, item.data.data(), item.size, item.blockCount, item.mode);
        if (specbleach1D_R)
            specbleach_load_noise_profile_for_mode(specbleach1D_R, item.data.data(), item.size, item.blockCount, item.mode);
    }

    for (const auto& item : profiles2D)
    {
        if (specbleach2D_L)
            specbleach_2d_load_noise_profile_for_mode(specbleach2D_L, item.data.data(), item.size, item.blockCount, item.mode);
        if (specbleach2D_R)
            specbleach_2d_load_noise_profile_for_mode(specbleach2D_R, item.data.data(), item.size, item.blockCount, item.mode);
    }

    // Synchronize loaded profiles between engines
    syncNoiseProfiles(0);
    syncNoiseProfiles(1);
}

void NoiseRepellentAudioProcessor::syncNoiseProfiles(int sourceAlgoMode)
{
    auto syncHandle = [](auto getProfileFn, auto getSizeFn, auto getBlockCountFn, auto isAvailableFn,
                         auto loadFn, auto targetHandle, auto sourceHandle) {
        if (!sourceHandle || !targetHandle) return;

        uint32_t sz = getSizeFn(sourceHandle);
        if (sz == 0) return;

        // Find a fallback profile (first available mode in source handle)
        float* fallbackP = nullptr;
        uint32_t fallbackBc = 0;

        for (int m = 1; m <= 4; ++m) {
            if (isAvailableFn(sourceHandle, m)) {
                fallbackP = getProfileFn(sourceHandle, m);
                fallbackBc = getBlockCountFn(sourceHandle, m);
                break;
            }
        }

        if (!fallbackP) return; // No profile learned in source handle yet

        // Load available mode profiles or use fallback so no target mode is left empty
        for (int mode = 1; mode <= 4; ++mode) {
            if (isAvailableFn(sourceHandle, mode)) {
                float* p = getProfileFn(sourceHandle, mode);
                uint32_t bc = getBlockCountFn(sourceHandle, mode);
                if (p) {
                    loadFn(targetHandle, p, sz, bc, mode);
                }
            } else {
                loadFn(targetHandle, fallbackP, sz, fallbackBc, mode);
            }
        }
    };

    if (sourceAlgoMode == 0)
    {
        // Sync 1D -> 2D
        syncHandle(specbleach_get_noise_profile_for_mode,
                   specbleach_get_noise_profile_size,
                   specbleach_get_noise_profile_block_count_for_mode,
                   specbleach_noise_profile_available_for_mode,
                   specbleach_2d_load_noise_profile_for_mode,
                   specbleach2D_L, specbleach1D_L);

        syncHandle(specbleach_get_noise_profile_for_mode,
                   specbleach_get_noise_profile_size,
                   specbleach_get_noise_profile_block_count_for_mode,
                   specbleach_noise_profile_available_for_mode,
                   specbleach_2d_load_noise_profile_for_mode,
                   specbleach2D_R, specbleach1D_R);
    }
    else if (sourceAlgoMode == 1)
    {
        // Sync 2D -> 1D
        syncHandle(specbleach_2d_get_noise_profile_for_mode,
                   specbleach_2d_get_noise_profile_size,
                   specbleach_2d_get_noise_profile_block_count_for_mode,
                   specbleach_2d_noise_profile_available_for_mode,
                   specbleach_load_noise_profile_for_mode,
                   specbleach1D_L, specbleach2D_L);

        syncHandle(specbleach_2d_get_noise_profile_for_mode,
                   specbleach_2d_get_noise_profile_size,
                   specbleach_2d_get_noise_profile_block_count_for_mode,
                   specbleach_2d_noise_profile_available_for_mode,
                   specbleach_load_noise_profile_for_mode,
                   specbleach1D_R, specbleach2D_R);
    }
}

void NoiseRepellentAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    ensureEnginesInitialized(sampleRate);

    // Set initial latency based on current algorithm mode
    currentAlgoMode = static_cast<int>(parameters.getRawParameterValue("algorithm_mode")->load());
    uint32_t latency = 0;
    if (currentAlgoMode == 0 && specbleach1D_L) {
        latency = specbleach_get_latency(specbleach1D_L);
    } else if (currentAlgoMode == 1 && specbleach2D_L) {
        latency = specbleach_2d_get_latency(specbleach2D_L);
    }
    setLatencySamples(static_cast<int>(latency));

    juce::dsp::ProcessSpec spec;
    spec.sampleRate = sampleRate;
    spec.maximumBlockSize = static_cast<juce::uint32>(samplesPerBlock);
    spec.numChannels = static_cast<juce::uint32>(getTotalNumOutputChannels());

    dryWetMixer.prepare(spec);
    dryWetMixer.setMixingRule(juce::dsp::DryWetMixingRule::linear);
    dryWetMixer.setWetLatency(static_cast<float>(latency));

    // Prepare buffer and state for engine crossfading
    crossfadeBuffer.setSize(getTotalNumOutputChannels(), samplesPerBlock, false, false, true);
    crossfadeProgress = 1.0f;
    targetAlgoMode = currentAlgoMode;

    // Reset FFT accumulation
    fftAccumInput.fill(0.0f);
    fftAccumOutput.fill(0.0f);
    fftAccumCount = 0;
}

void NoiseRepellentAudioProcessor::releaseResources()
{
    dryWetMixer.reset();
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

    const float smoothing = parameters.getRawParameterValue("smoothing_factor")->load();
    const float nlmSmoothing = parameters.getRawParameterValue("nlm_smoothing")->load();

    // Masking: apply cubic curve mapping (matches original LV2 behavior)
    const float maskingRaw = parameters.getRawParameterValue("masking_depth")->load() / 100.0f;
    const float masking = 1.0f - std::pow(1.0f - maskingRaw, 3.0f);

    const float whitening = parameters.getRawParameterValue("whitening_factor")->load();

    // Suppression: pass raw 0-100 value — libspecbleach divides by 100 internally
    const float suppression = parameters.getRawParameterValue("suppression_strength")->load();

    const bool residualListen = parameters.getRawParameterValue("residual_listen")->load() > 0.5f;

    // Sync profiles if manual noise learning was just turned off
    if (wasLearning && !learnNoise) {
        syncNoiseProfiles(currentAlgoMode);
    }
    wasLearning = learnNoise;

    // Update latency dynamically, sync profiles, and start smooth crossfade when algorithm mode changes
    if (algoMode != currentAlgoMode) {
        syncNoiseProfiles(currentAlgoMode);
        targetAlgoMode = algoMode;
        currentAlgoMode = algoMode;

        // Initiate 30 ms smooth crossfade transition to prevent clicks/pops
        int crossfadeSamples = static_cast<int>(currentSampleRate * 0.030);
        if (crossfadeSamples < 1) crossfadeSamples = 1;
        crossfadeProgress = 0.0f;
        crossfadeStep = 1.0f / static_cast<float>(crossfadeSamples);

        uint32_t latency = 0;
        if (algoMode == 0 && specbleach1D_L) {
            latency = specbleach_get_latency(specbleach1D_L);
        } else if (algoMode == 1 && specbleach2D_L) {
            latency = specbleach_2d_get_latency(specbleach2D_L);
        }
        setLatencySamples(static_cast<int>(latency));
        dryWetMixer.setWetLatency(static_cast<float>(latency));
    }

    // Save dry input copy for FFT visualization before processing
    std::vector<float> dryInputL(static_cast<size_t>(numSamples));
    if (numChannels >= 1)
        std::copy_n(buffer.getReadPointer(0), numSamples, dryInputL.begin());

    // Push dry samples into JUCE DryWetMixer before in-place DSP processing
    juce::dsp::AudioBlock<float> audioBlock(buffer);
    dryWetMixer.pushDrySamples(audioBlock);

    // Prepare parameters for DSP engines
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
    p.tonal_reduction = tonalRed;

    SpectralBleach2DDenoiserParameters p2;
    p2.learn_noise = learnNoise ? 1 : 0;
    p2.residual_listen = residualListen;
    p2.reduction_amount = masterRed;
    p2.smoothing_factor = nlmSmoothing;
    p2.whitening_factor = whitening;
    p2.adaptive_noise = adaptiveNoise ? 1 : 0;
    p2.noise_estimation_method = adaptiveMethod;
    p2.nlm_masking_protection = masking;
    p2.suppression_strength = suppression;
    p2.aggressiveness = aggressiveness;
    p2.tonal_reduction = tonalRed;

    if (crossfadeProgress < 1.0f)
    {
        // Smooth crossfade mode: run both engines and blend sample-by-sample
        if (crossfadeBuffer.getNumChannels() < numChannels || crossfadeBuffer.getNumSamples() < numSamples)
            crossfadeBuffer.setSize(numChannels, numSamples, false, false, true);

        for (int ch = 0; ch < numChannels; ++ch)
            std::copy_n(buffer.getReadPointer(ch), numSamples, crossfadeBuffer.getWritePointer(ch));

        // Process 1D on buffer
        if (specbleach1D_L) {
            specbleach_load_parameters(specbleach1D_L, p);
            float* ch0 = buffer.getWritePointer(0);
            specbleach_process(specbleach1D_L, static_cast<uint32_t>(numSamples), ch0, ch0);
        }
        if (numChannels >= 2 && specbleach1D_R) {
            specbleach_load_parameters(specbleach1D_R, p);
            float* ch1 = buffer.getWritePointer(1);
            specbleach_process(specbleach1D_R, static_cast<uint32_t>(numSamples), ch1, ch1);
        }

        // Process 2D on crossfadeBuffer
        if (specbleach2D_L) {
            specbleach_2d_load_parameters(specbleach2D_L, p2);
            float* ch0 = crossfadeBuffer.getWritePointer(0);
            specbleach_2d_process(specbleach2D_L, static_cast<uint32_t>(numSamples), ch0, ch0);
        }
        if (numChannels >= 2 && specbleach2D_R) {
            specbleach_2d_load_parameters(specbleach2D_R, p2);
            float* ch1 = crossfadeBuffer.getWritePointer(1);
            specbleach_2d_process(specbleach2D_R, static_cast<uint32_t>(numSamples), ch1, ch1);
        }

        // Equal-power crossfade blend into buffer
        for (int ch = 0; ch < numChannels; ++ch) {
            float* out1D = buffer.getWritePointer(ch);
            const float* out2D = crossfadeBuffer.getReadPointer(ch);

            float currentP = crossfadeProgress;
            for (int s = 0; s < numSamples; ++s) {
                float pVal = std::min(1.0f, currentP);
                float rad = pVal * juce::MathConstants<float>::halfPi;
                float w1D = (targetAlgoMode == 1) ? std::cos(rad) : std::sin(rad);
                float w2D = (targetAlgoMode == 1) ? std::sin(rad) : std::cos(rad);

                out1D[s] = w1D * out1D[s] + w2D * out2D[s];
                currentP += crossfadeStep;
            }
        }
        crossfadeProgress += numSamples * crossfadeStep;
        if (crossfadeProgress > 1.0f) crossfadeProgress = 1.0f;
    }
    else if (algoMode == 0)
    {
        // Steady-state 1D mode
        if (specbleach1D_L) {
            specbleach_load_parameters(specbleach1D_L, p);
            float* ch0 = buffer.getWritePointer(0);
            specbleach_process(specbleach1D_L, static_cast<uint32_t>(numSamples), ch0, ch0);
        }
        if (numChannels >= 2 && specbleach1D_R) {
            specbleach_load_parameters(specbleach1D_R, p);
            float* ch1 = buffer.getWritePointer(1);
            specbleach_process(specbleach1D_R, static_cast<uint32_t>(numSamples), ch1, ch1);
        }
    }
    else if (algoMode == 1)
    {
        // Steady-state 2D mode
        if (specbleach2D_L) {
            specbleach_2d_load_parameters(specbleach2D_L, p2);
            float* ch0 = buffer.getWritePointer(0);
            specbleach_2d_process(specbleach2D_L, static_cast<uint32_t>(numSamples), ch0, ch0);
        }
        if (numChannels >= 2 && specbleach2D_R) {
            specbleach_2d_load_parameters(specbleach2D_R, p2);
            float* ch1 = buffer.getWritePointer(1);
            specbleach_2d_process(specbleach2D_R, static_cast<uint32_t>(numSamples), ch1, ch1);
        }
    }

    // Apply Soft Crossfade Bypass using JUCE DryWetMixer (with dry latency compensation)
    dryWetMixer.setWetMixProportion(isBypassed ? 0.0f : 1.0f);
    dryWetMixer.mixWetSamples(audioBlock);

    // ── FFT-based spectral frame for GUI visualization ──
    // Accumulate samples until we have a full FFT window (kFftSize samples)
    const float* inputSrc = dryInputL.data();
    const float* outputSrc = (numChannels >= 1) ? buffer.getReadPointer(0) : nullptr;

    for (int s = 0; s < numSamples; ++s)
    {
        if (fftAccumCount < kFftSize) {
            fftAccumInput[fftAccumCount] = inputSrc[s];
            fftAccumOutput[fftAccumCount] = outputSrc ? outputSrc[s] : 0.0f;
            fftAccumCount++;
        }

        if (fftAccumCount >= kFftSize)
        {
            // We have a full FFT frame — compute and push to ring buffer
            int start1, size1, start2, size2;
            spectralFifo.prepareToWrite(1, start1, size1, start2, size2);

            if (size1 > 0)
            {
                SpectralFrame& frame = spectralBuffer[static_cast<size_t>(start1)];

                // ── Input FFT ──
                std::memcpy(fftInputWork.data(), fftAccumInput.data(), kFftSize * sizeof(float));
                std::fill(fftInputWork.begin() + kFftSize, fftInputWork.end(), 0.0f);
                fftWindow.multiplyWithWindowingTable(fftInputWork.data(), kFftSize);
                fftAnalyzer.performFrequencyOnlyForwardTransform(fftInputWork.data());

                for (size_t i = 0; i < kFftBins; ++i) {
                    float mag = fftInputWork[i] / static_cast<float>(kFftBins);
                    frame.inputMagnitudeDB[i] = 20.0f * std::log10(std::max(mag, 1e-7f));
                }

                // ── Output FFT ──
                std::memcpy(fftOutputWork.data(), fftAccumOutput.data(), kFftSize * sizeof(float));
                std::fill(fftOutputWork.begin() + kFftSize, fftOutputWork.end(), 0.0f);
                fftWindow.multiplyWithWindowingTable(fftOutputWork.data(), kFftSize);
                fftAnalyzer.performFrequencyOnlyForwardTransform(fftOutputWork.data());

                for (size_t i = 0; i < kFftBins; ++i) {
                    float mag = fftOutputWork[i] / static_cast<float>(kFftBins);
                    frame.outputMagnitudeDB[i] = 20.0f * std::log10(std::max(mag, 1e-7f));
                }

                // ── Noise profile from libspecbleach (stationary or live adaptive) ──
                const float* actualNoiseProfile = nullptr;
                uint32_t profileSize = 0;
                bool profileAvailable = false;

                bool isAdaptive = adaptiveNoise;

                if (algoMode == 0 && specbleach1D_L &&
                    (isAdaptive || specbleach_noise_profile_available_for_mode(specbleach1D_L, 1))) {
                    actualNoiseProfile = specbleach_get_active_noise_profile(specbleach1D_L);
                    profileSize = specbleach_get_noise_profile_size(specbleach1D_L);
                    profileAvailable = (actualNoiseProfile != nullptr && profileSize > 0);
                } else if (algoMode == 1 && specbleach2D_L &&
                           (isAdaptive || specbleach_2d_noise_profile_available_for_mode(specbleach2D_L, 1))) {
                    actualNoiseProfile = specbleach_2d_get_active_noise_profile(specbleach2D_L);
                    profileSize = specbleach_2d_get_noise_profile_size(specbleach2D_L);
                    profileAvailable = (actualNoiseProfile != nullptr && profileSize > 0);
                }

                frame.hasNoiseProfile = profileAvailable;
                frame.isLinked = (parameters.getRawParameterValue("link_reduction")->load() > 0.5f);
                frame.tonalPeaksHz.clear();

                if (profileAvailable && actualNoiseProfile) {
                    // In 2D mode, noise_profile_size is SB_PROCESSOR_CORE_FULL_FFT_SPECTRUM (e.g. 3008 bins).
                    // The real unique spectrum from DC (0 Hz) to Nyquist (Fs/2) is the first (profileSize / 2) + 1 bins.
                    size_t realProfileBins = profileSize;

                    const float maxProfileIdx = static_cast<float>(realProfileBins - 1);
                    const float maxFftIdx = static_cast<float>(kFftBins - 1);
                    // FFTW unnormalized power scaling offset: 20 * log10(N/2)
                    const float dbOffset = (maxProfileIdx > 0.0f) ? (20.0f * std::log10(maxProfileIdx)) : 0.0f;

                    for (size_t i = 0; i < kFftBins; ++i) {
                        float normPos = static_cast<float>(i) / maxFftIdx;
                        float exactP = normPos * maxProfileIdx;

                        size_t p0 = std::clamp(static_cast<size_t>(exactP), static_cast<size_t>(0), realProfileBins - 1);
                        size_t p1 = std::min(p0 + 1, realProfileBins - 1);
                        float frac = exactP - static_cast<float>(p0);

                        float val0 = actualNoiseProfile[p0];
                        float val1 = actualNoiseProfile[p1];
                        float interpVal = (1.0f - frac) * val0 + frac * val1;

                        float rawDb = 10.0f * std::log10(std::max(interpVal, 1e-12f));
                        frame.noiseFloorDB[i] = rawDb - dbOffset;
                    }

                    // Query libspecbleach directly for reported tonal peak frequencies computed on this real noise profile
                    std::array<float, 16> peakBuf{};
                    uint32_t numPeaks = 0;
                    if (algoMode == 0 && specbleach1D_L) {
                        numPeaks = specbleach_get_tonal_peaks_for_profile(specbleach1D_L, actualNoiseProfile, static_cast<uint32_t>(realProfileBins), peakBuf.data(), static_cast<uint32_t>(peakBuf.size()));
                    } else if (algoMode == 1 && specbleach2D_L) {
                        numPeaks = specbleach_2d_get_tonal_peaks_for_profile(specbleach2D_L, actualNoiseProfile, static_cast<uint32_t>(realProfileBins), peakBuf.data(), static_cast<uint32_t>(peakBuf.size()));
                    }

                    for (uint32_t i = 0; i < numPeaks; ++i) {
                        frame.tonalPeaksHz.push_back(peakBuf[i]);
                    }
                } else {
                    frame.noiseFloorDB.fill(-120.0f);
                }

                spectralFifo.finishedWrite(1);
            }

            // Shift: keep last quarter for overlap (hop = 75%)
            constexpr size_t kHop = kFftSize / 4;
            std::memmove(fftAccumInput.data(), fftAccumInput.data() + kHop, (kFftSize - kHop) * sizeof(float));
            std::memmove(fftAccumOutput.data(), fftAccumOutput.data() + kHop, (kFftSize - kHop) * sizeof(float));
            fftAccumCount = kFftSize - kHop;
        }
    }
}

void NoiseRepellentAudioProcessor::processBlockBypassed(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& /*midiMessages*/)
{
    const int numSamples = buffer.getNumSamples();
    const int numChannels = buffer.getNumChannels();

    if (numSamples == 0 || numChannels == 0)
        return;

    juce::dsp::AudioBlock<float> audioBlock(buffer);
    dryWetMixer.pushDrySamples(audioBlock);
    dryWetMixer.setWetMixProportion(0.0f);
    dryWetMixer.mixWetSamples(audioBlock);
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
    if (specbleach1D_L) specbleach_reset_noise_profile(specbleach1D_L);
    if (specbleach1D_R) specbleach_reset_noise_profile(specbleach1D_R);
    if (specbleach2D_L) specbleach_2d_reset_noise_profile(specbleach2D_L);
    if (specbleach2D_R) specbleach_2d_reset_noise_profile(specbleach2D_R);
}

bool NoiseRepellentAudioProcessor::hasNoiseProfile() const
{
    auto* paramAlgo = parameters.getRawParameterValue("algorithm_mode");
    if (!paramAlgo)
        return false;

    int algoMode = static_cast<int>(paramAlgo->load());

    if (algoMode == 0 && specbleach1D_L)
    {
        for (int mode = 1; mode <= 4; ++mode)
        {
            if (specbleach_noise_profile_available_for_mode(specbleach1D_L, mode))
                return true;
        }
    }
    else if (algoMode == 1 && specbleach2D_L)
    {
        for (int mode = 1; mode <= 4; ++mode)
        {
            if (specbleach_2d_noise_profile_available_for_mode(specbleach2D_L, mode))
                return true;
        }
    }

    return false;
}

juce::AudioProcessorEditor* NoiseRepellentAudioProcessor::createEditor()
{
    return new NoiseRepellentAudioProcessorEditor(*this);
}

void NoiseRepellentAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
    auto state = parameters.copyState();

    juce::ValueTree profilesTree("LEARNED_PROFILES");

    auto saveProfiles = [&](SpectralBleachHandle instance, const char* tagName,
                            auto getProfileFn, auto getSizeFn, auto getBlockCountFn, auto isAvailableFn) {
        // Modes 1-4: ROLLING_MEAN=1, MEDIAN=2, MAX=3, MINIMUM=4
        for (int mode = 1; mode <= 4; ++mode) {
            if (instance && isAvailableFn(instance, mode)) {
                float* profile = getProfileFn(instance, mode);
                uint32_t profileSize = getSizeFn(instance);
                uint32_t blockCount = getBlockCountFn(instance, mode);

                if (profile && profileSize > 0) {
                    juce::MemoryBlock mb(profile, profileSize * sizeof(float));
                    juce::ValueTree node(tagName);
                    node.setProperty("mode", mode, nullptr);
                    node.setProperty("size", static_cast<int>(profileSize), nullptr);
                    node.setProperty("blockCount", static_cast<int>(blockCount), nullptr);
                    node.setProperty("data", mb.toBase64Encoding(), nullptr);
                    profilesTree.appendChild(node, nullptr);
                }
            }
        }
    };

    saveProfiles(specbleach1D_L, "PROFILE_1D",
                 specbleach_get_noise_profile_for_mode,
                 specbleach_get_noise_profile_size,
                 specbleach_get_noise_profile_block_count_for_mode,
                 specbleach_noise_profile_available_for_mode);

    saveProfiles(specbleach2D_L, "PROFILE_2D",
                 specbleach_2d_get_noise_profile_for_mode,
                 specbleach_2d_get_noise_profile_size,
                 specbleach_2d_get_noise_profile_block_count_for_mode,
                 specbleach_2d_noise_profile_available_for_mode);

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

                if (size > 0 && base64Data.isNotEmpty() && mode >= 1 && mode <= 4)
                {
                    juce::MemoryBlock mb;
                    if (mb.fromBase64Encoding(base64Data) && mb.getSize() >= size * sizeof(float))
                    {
                        const float* floatArray = reinterpret_cast<const float*>(mb.getData());
                        if (node.hasType("PROFILE_1D"))
                        {
                            if (specbleach1D_L)
                                specbleach_load_noise_profile_for_mode(specbleach1D_L, floatArray, size, blockCount, mode);
                            if (specbleach1D_R)
                                specbleach_load_noise_profile_for_mode(specbleach1D_R, floatArray, size, blockCount, mode);
                        }
                        else if (node.hasType("PROFILE_2D"))
                        {
                            if (specbleach2D_L)
                                specbleach_2d_load_noise_profile_for_mode(specbleach2D_L, floatArray, size, blockCount, mode);
                            if (specbleach2D_R)
                                specbleach_2d_load_noise_profile_for_mode(specbleach2D_R, floatArray, size, blockCount, mode);
                        }
                    }
                }
            }

            // Sync loaded profiles to ensure both engines are fully populated
            syncNoiseProfiles(0);
            syncNoiseProfiles(1);
        }
    }
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new NoiseRepellentAudioProcessor();
}
