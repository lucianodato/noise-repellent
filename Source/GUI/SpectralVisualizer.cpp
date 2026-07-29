/*
noise-repellent -- Noise Reduction JUCE Plugin

Copyright 2026 Luciano Dato <lucianodato@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.
*/

#include "SpectralVisualizer.h"
#include "LookAndFeel.h"

SpectralVisualizerComponent::SpectralVisualizerComponent(NoiseRepellentAudioProcessor& p)
    : processor(p)
{
    startTimerHz(60); // 60 FPS visualizer
}

SpectralVisualizerComponent::~SpectralVisualizerComponent()
{
    stopTimer();
}

void SpectralVisualizerComponent::timerCallback()
{
    if (processor.getNextSpectralFrame(currentFrame))
    {
        repaint();
    }
}

void SpectralVisualizerComponent::resized()
{
}

void SpectralVisualizerComponent::paint(juce::Graphics& g)
{
    g.fillAll(juce::Colour(0xff090a0c));

    const float w = static_cast<float>(getWidth());
    const float h = static_cast<float>(getHeight());

    // Grid lines
    g.setColour(juce::Colour(0xff181c26));
    g.setFont(10.0f);

    // dB Y-grid (-80 dB to +10 dB)
    for (int db = -80; db <= 10; db += 20)
    {
        float y = h - ((static_cast<float>(db) + 90.0f) / 100.0f) * h;
        g.drawHorizontalLine(static_cast<int>(y), 0.0f, w);

        g.setColour(juce::Colour(0xff505a6e));
        juce::String label = (db > 0 ? "+" : "") + juce::String(db) + " dB";
        g.drawText(label, 8, static_cast<int>(y) - 12, 60, 12, juce::Justification::left);
        g.setColour(juce::Colour(0xff181c26));
    }

    // Frequency X-grid (50Hz - 20kHz logarithmic)
    static const float freqs[] = { 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000 };
    for (float f : freqs)
    {
        float x = (std::log10(f / 20.0f) / std::log10(20000.0f / 20.0f)) * w;
        g.drawVerticalLine(static_cast<int>(x), 0.0f, h);

        g.setColour(juce::Colour(0xff505a6e));
        juce::String label = f >= 1000.0f ? juce::String(f / 1000.0f, 0) + "k" : juce::String(static_cast<int>(f));
        g.drawText(label, static_cast<int>(x) + 4, static_cast<int>(h) - 14, 40, 12, juce::Justification::left);
        g.setColour(juce::Colour(0xff181c26));
    }

    // 1. Tonal Peak Markers (Muted Violet) - Only draw if detected
    if (!currentFrame.tonalPeaksHz.empty())
    {
        for (float peakHz : currentFrame.tonalPeaksHz)
        {
            float x = (std::log10(peakHz / 20.0f) / std::log10(20000.0f / 20.0f)) * w;

            juce::Line<float> line(x, 0.0f, x, h);
            float dashLengths[] = { 3.0f, 3.0f };
            g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks);
            g.drawDashedLine(line, dashLengths, 2, 1.5f);

            g.fillRoundedRectangle(x - 16.0f, 8.0f, 32.0f, 14.0f, 3.0f);
            g.setColour(juce::Colours::white);
            g.setFont(9.0f);
            juce::String tag = peakHz >= 1000.0f ? juce::String(peakHz / 1000.0f, 1) + "kHz" : juce::String(static_cast<int>(peakHz)) + "Hz";
            g.drawText(tag, static_cast<int>(x) - 16, 8, 32, 14, juce::Justification::centred);
        }
    }

    // 2. Noise Floor Profile Curve (Solid Warm Amber) - Only draw if noise has been captured
    if (currentFrame.hasNoiseProfile)
    {
        juce::Path noisePath;
        for (size_t i = 0; i < NoiseRepellentAudioProcessor::kFftSize; ++i)
        {
            float normX = static_cast<float>(i) / static_cast<float>(NoiseRepellentAudioProcessor::kFftSize);
            float x = normX * w;
            float level = currentFrame.noiseFloor[i];
            float y = h * (1.0f - level);

            if (i == 0) noisePath.startNewSubPath(x, y);
            else noisePath.lineTo(x, y);
        }
        g.setColour(NoiseRepellentLookAndFeel::kColorNoiseProfile);
        g.strokePath(noisePath, juce::PathStrokeType(2.0f));
    }

    // 3. Input Signal Curve (Soft Muted Slate)
    juce::Path inputPath;
    for (size_t i = 0; i < NoiseRepellentAudioProcessor::kFftSize; ++i)
    {
        float normX = static_cast<float>(i) / static_cast<float>(NoiseRepellentAudioProcessor::kFftSize);
        float x = normX * w;
        float level = std::min(1.0f, currentFrame.inputMagnitude[i]);
        float y = h * (1.0f - level);

        if (i == 0) inputPath.startNewSubPath(x, y);
        else inputPath.lineTo(x, y);
    }
    g.setColour(NoiseRepellentLookAndFeel::kColorInputSignal);
    g.strokePath(inputPath, juce::PathStrokeType(1.5f));

    // 4. Denoised Output Area (Solid Translucent Slate Blue)
    juce::Path outputPath;
    outputPath.startNewSubPath(0.0f, h);
    for (size_t i = 0; i < NoiseRepellentAudioProcessor::kFftSize; ++i)
    {
        float normX = static_cast<float>(i) / static_cast<float>(NoiseRepellentAudioProcessor::kFftSize);
        float x = normX * w;
        float level = std::min(1.0f, currentFrame.outputMagnitude[i]);
        float y = h * (1.0f - level);
        outputPath.lineTo(x, y);
    }
    outputPath.lineTo(w, h);
    outputPath.closeSubPath();

    g.setColour(NoiseRepellentLookAndFeel::kColorDenoising.withAlpha(0.25f));
    g.fillPath(outputPath);
}
