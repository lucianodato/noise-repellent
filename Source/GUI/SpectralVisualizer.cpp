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

#include "SpectralVisualizer.h"
#include "LookAndFeel.h"
#include <cmath>
#include <algorithm>

SpectralVisualizerComponent::SpectralVisualizerComponent(NoiseRepellentAudioProcessor& p)
    : processor(p)
{
    smoothedInputDB.fill(-100.0f);
    smoothedOutputDB.fill(-100.0f);
    startTimerHz(30);
}

SpectralVisualizerComponent::~SpectralVisualizerComponent()
{
    stopTimer();
}

void SpectralVisualizerComponent::timerCallback()
{
    if (processor.getNextSpectralFrame(currentFrame))
    {
        const size_t numBins = NoiseRepellentAudioProcessor::kFftBins;

        if (!isSmoothedInitialized)
        {
            smoothedInputDB = currentFrame.inputMagnitudeDB;
            smoothedOutputDB = currentFrame.outputMagnitudeDB;
            isSmoothedInitialized = true;
        }
        else
        {
            // Asymmetric Exponential Moving Average (Fast Attack, Smooth Fluid Decay)
            constexpr float kAttackAlpha  = 0.35f; // Fast response to audio transients
            constexpr float kReleaseAlpha = 0.10f; // Smooth, liquid decay on falling signals

            for (size_t i = 0; i < numBins; ++i)
            {
                float targetIn = currentFrame.inputMagnitudeDB[i];
                float alphaIn = (targetIn > smoothedInputDB[i]) ? kAttackAlpha : kReleaseAlpha;
                smoothedInputDB[i] += alphaIn * (targetIn - smoothedInputDB[i]);

                float targetOut = currentFrame.outputMagnitudeDB[i];
                float alphaOut = (targetOut > smoothedOutputDB[i]) ? kAttackAlpha : kReleaseAlpha;
                smoothedOutputDB[i] += alphaOut * (targetOut - smoothedOutputDB[i]);
            }
        }

        repaint();
    }
}

void SpectralVisualizerComponent::resized()
{
}

// Map a frequency (Hz) to an X pixel position using logarithmic scale
static float freqToX(float freqHz, float width, float minFreq = 20.0f, float maxFreq = 20000.0f)
{
    if (freqHz <= minFreq) return 0.0f;
    if (freqHz >= maxFreq) return width;
    return (std::log10(freqHz / minFreq) / std::log10(maxFreq / minFreq)) * width;
}

// Map a dB value to a Y pixel position
static float dbToY(float db, float height, float minDB = -100.0f, float maxDB = -20.0f)
{
    float clamped = juce::jlimit(minDB, maxDB, db);
    float normalized = (clamped - minDB) / (maxDB - minDB); // 0 at bottom, 1 at top
    return height * (1.0f - normalized);
}

// Apply variable fractional-octave spatial frequency smoothing across FFT bins (UI display only).
// Low frequencies (bass/hum) get minimal smoothing (preserving sharp peaks),
// while high frequencies get progressively wider smoothing windows matching log-scale visuals.
static void smoothBinsLogarithmicFrequencyDomain(const float* src, float* dst, size_t numBins)
{
    if (numBins == 0) return;

    for (size_t i = 0; i < numBins; ++i)
    {
        // Half-window radius grows proportionally with bin index i (fractional octave)
        // i = 5 (~50 Hz)   -> radius = 0 (exact bin value)
        // i = 20 (~200 Hz) -> radius = 1 (3-bin window)
        // i = 100 (~1 kHz) -> radius = 3 (7-bin window)
        // i = 500 (~5 kHz) -> radius = 12 (25-bin window)
        // i = 1000 (~10 kHz)-> radius = 20 (41-bin window)
        int radius = std::clamp(static_cast<int>(std::round(static_cast<float>(i) * 0.025f)), 0, 20);

        if (radius == 0)
        {
            dst[i] = src[i];
            continue;
        }

        int minIdx = std::max(0, static_cast<int>(i) - radius);
        int maxIdx = std::min(static_cast<int>(numBins - 1), static_cast<int>(i) + radius);

        float sum = 0.0f;
        float totalWeight = 0.0f;

        for (int k = minIdx; k <= maxIdx; ++k)
        {
            // Triangular window weight: highest weight at center bin i, tapering down to edges
            float dist = std::abs(static_cast<float>(k - static_cast<int>(i)));
            float weight = 1.0f - (dist / static_cast<float>(radius + 1));

            sum += src[k] * weight;
            totalWeight += weight;
        }

        dst[i] = (totalWeight > 0.0f) ? (sum / totalWeight) : src[i];
    }
}

void SpectralVisualizerComponent::paint(juce::Graphics& g)
{
    g.fillAll(juce::Colour(0xff232832));

    const float w = static_cast<float>(getWidth());
    const float h = static_cast<float>(getHeight());
    const float minDB = -100.0f;
    const float maxDB = -20.0f;
    const float minFreq = 20.0f;
    const float maxFreq = 20000.0f;

    const double sampleRate = processor.getSampleRate();
    const size_t numBins = NoiseRepellentAudioProcessor::kFftBins;
    const float binWidth = static_cast<float>(sampleRate) / static_cast<float>(NoiseRepellentAudioProcessor::kFftSize);

    // Frequency-domain spatially smoothed display buffers (UI display smoothing only)
    std::vector<float> freqSmoothedInput(numBins);
    std::vector<float> freqSmoothedOutput(numBins);
    std::vector<float> freqSmoothedProfile(numBins);

    smoothBinsLogarithmicFrequencyDomain(smoothedInputDB.data(), freqSmoothedInput.data(), numBins);
    smoothBinsLogarithmicFrequencyDomain(smoothedOutputDB.data(), freqSmoothedOutput.data(), numBins);
    if (currentFrame.hasNoiseProfile) {
        smoothBinsLogarithmicFrequencyDomain(currentFrame.noiseFloorDB.data(), freqSmoothedProfile.data(), numBins);
    }

    // ── Grid lines ──
    g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));

    // dB Y-grid (-90 dB to -20 dB in 10 dB steps)
    for (int db = -90; db <= -20; db += 10)
    {
        float y = dbToY(static_cast<float>(db), h, minDB, maxDB);
        g.setColour(juce::Colour(0xff3d4657));
        g.drawHorizontalLine(static_cast<int>(y), 0.0f, w);

        g.setColour(juce::Colour(0xffa8b3c4));
        juce::String label = juce::String(db) + " dB";
        g.drawText(label, 8, static_cast<int>(y) - 12, 60, 12, juce::Justification::left);
    }

    // Frequency X-grid (logarithmic)
    static const float freqs[] = { 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000 };
    for (float f : freqs)
    {
        float x = freqToX(f, w, minFreq, maxFreq);
        g.setColour(juce::Colour(0xff3d4657));
        g.drawVerticalLine(static_cast<int>(x), 0.0f, h);

        g.setColour(juce::Colour(0xffa8b3c4));
        juce::String label = f >= 1000.0f ? juce::String(f / 1000.0f, 0) + "k" : juce::String(static_cast<int>(f));
        g.drawText(label, static_cast<int>(x) + 4, static_cast<int>(h) - 14, 40, 12, juce::Justification::left);
    }

    // ── Helper: build a path by mapping FFT bins to log-frequency X positions ──
    auto buildLogFreqPath = [&](const float* dbData, size_t bins) -> juce::Path {
        juce::Path path;
        bool started = false;

        for (size_t i = 1; i < bins; ++i) // skip bin 0 (DC)
        {
            float freq = static_cast<float>(i) * binWidth;
            if (freq < minFreq || freq > maxFreq) continue;

            float x = freqToX(freq, w, minFreq, maxFreq);
            float y = dbToY(dbData[i], h, minDB, maxDB);

            if (!started) {
                path.startNewSubPath(x, y);
                started = true;
            } else {
                path.lineTo(x, y);
            }
        }
        return path;
    };

    // 1. Input Signal Spectrum (Filled Translucent Area)
    {
        juce::Path inputAreaPath;
        bool started = false;
        float firstX = 0.0f;
        float lastX = 0.0f;

        for (size_t i = 1; i < numBins; ++i)
        {
            float freq = static_cast<float>(i) * binWidth;
            if (freq < minFreq || freq > maxFreq) continue;

            float x = freqToX(freq, w, minFreq, maxFreq);
            float y = dbToY(freqSmoothedInput[i], h, minDB, maxDB);

            if (!started) {
                firstX = x;
                inputAreaPath.startNewSubPath(x, h);
                inputAreaPath.lineTo(x, y);
                started = true;
            } else {
                inputAreaPath.lineTo(x, y);
            }
            lastX = x;
        }

        if (started) {
            inputAreaPath.lineTo(lastX, h);
            inputAreaPath.lineTo(firstX, h);
            inputAreaPath.closeSubPath();

            g.setColour(NoiseRepellentLookAndFeel::kColorInputSignal.withAlpha(0.30f));
            g.fillPath(inputAreaPath);
        }
    }

    // 2. Denoised Output Signal Curve (Bright Solid Cyan Line)
    {
        juce::Path outputPath = buildLogFreqPath(freqSmoothedOutput.data(), numBins);
        g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
        g.strokePath(outputPath, juce::PathStrokeType(1.8f));
    }

    // 3. Noise Floor Profile Curve (Solid Warm Amber Line)
    if (currentFrame.hasNoiseProfile)
    {
        juce::Path noisePath = buildLogFreqPath(freqSmoothedProfile.data(), numBins);
        g.setColour(NoiseRepellentLookAndFeel::kColorNoiseProfile);
        g.strokePath(noisePath, juce::PathStrokeType(2.0f));
    }

    // 4. Tonal Peak Markers (Dashed vertical lines with staggered multi-tier frequency tags; shown ONLY when Reduction is UNLINKED)
    if (!currentFrame.isLinked && currentFrame.hasNoiseProfile && !currentFrame.tonalPeaksHz.empty())
    {
        // 3 Y-tiers starting at Y = 42.0f to guarantee zero collision with top HUD overlay bar
        constexpr float kStartTagY = 42.0f;
        float lastTagRight[3] = { -100.0f, -100.0f, -100.0f };

        for (float peakHz : currentFrame.tonalPeaksHz)
        {
            float x = freqToX(peakHz, w, minFreq, maxFreq);
            if (x <= 4.0f || x >= w - 4.0f) continue;

            // Dashed vertical line across the display
            juce::Line<float> line(x, 0.0f, x, h);
            float dashLengths[] = { 3.0f, 3.0f };
            g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks.withAlpha(0.75f));
            g.drawDashedLine(line, dashLengths, 2, 1.5f);

            // Frequency Tag Badge
            juce::String tag = peakHz >= 1000.0f ? juce::String(peakHz / 1000.0f, 1) + "k" : juce::String(static_cast<int>(peakHz)) + "Hz";
            float tagW = 38.0f;
            float tagH = 14.0f;
            float tagX = std::clamp(x - tagW * 0.5f, 4.0f, w - tagW - 4.0f);

            // Find first non-colliding Y-tier
            int chosenTier = -1;
            for (int t = 0; t < 3; ++t) {
                if (tagX > lastTagRight[t] + 2.0f) {
                    chosenTier = t;
                    break;
                }
            }

            if (chosenTier != -1)
            {
                float tagY = kStartTagY + static_cast<float>(chosenTier) * 16.0f;

                // Draw vertical leader connector line from top to badge
                g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks.withAlpha(0.5f));
                g.drawVerticalLine(static_cast<int>(x), 0.0f, tagY);

                g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks.withAlpha(0.9f));
                g.fillRoundedRectangle(tagX, tagY, tagW, tagH, 3.0f);

                g.setColour(juce::Colours::white);
                g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));
                g.drawText(tag, static_cast<int>(tagX), static_cast<int>(tagY), static_cast<int>(tagW), static_cast<int>(tagH), juce::Justification::centred);

                lastTagRight[chosenTier] = tagX + tagW;
            }
        }
    }

    // 5. HUD Color Legend (Top-Right Overlay, positioned to the left of Advanced Controls button)
    {
        const bool showTonalSwatch = !currentFrame.isLinked;
        const float legendW = showTonalSwatch ? 292.0f : 205.0f;
        const float legendH = 24.0f;
        const float advButtonW = 148.0f;
        const float legendX = std::max(250.0f, w - advButtonW - 12.0f - 10.0f - legendW);
        const float legendY = 10.0f;

        // Semi-transparent dark glass background
        g.setColour(juce::Colour(0xeb252a35));
        g.fillRoundedRectangle(legendX, legendY, legendW, legendH, 4.0f);
        g.setColour(juce::Colour(0xff4c566a));
        g.drawRoundedRectangle(legendX, legendY, legendW, legendH, 4.0f, 1.0f);

        g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel, juce::Font::bold));

        float curX = legendX + 10.0f;

        // Swatch 1: Live Input (Filled Area)
        g.setColour(NoiseRepellentLookAndFeel::kColorInputSignal.withAlpha(0.70f));
        g.fillRect(curX, legendY + 7.0f, 11.0f, 8.0f);
        curX += 15.0f;
        g.setColour(juce::Colour(0xffd8e0ec));
        g.drawText("Input", static_cast<int>(curX), static_cast<int>(legendY), 32, static_cast<int>(legendH), juce::Justification::left);
        curX += 46.0f;

        // Swatch 2: Noise Profile (Solid Amber Line)
        g.setColour(NoiseRepellentLookAndFeel::kColorNoiseProfile);
        g.drawLine(curX, legendY + 11.0f, curX + 13.0f, legendY + 11.0f, 2.0f);
        curX += 17.0f;
        g.setColour(juce::Colour(0xffd8e0ec));
        g.drawText("Profile", static_cast<int>(curX), static_cast<int>(legendY), 36, static_cast<int>(legendH), juce::Justification::left);
        curX += 50.0f;

        // Swatch 3: Output Area (Solid Cyan Line)
        g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
        g.drawLine(curX, legendY + 11.0f, curX + 13.0f, legendY + 11.0f, 2.0f);
        curX += 17.0f;
        g.setColour(juce::Colour(0xffd8e0ec));
        g.drawText("Output", static_cast<int>(curX), static_cast<int>(legendY), 38, static_cast<int>(legendH), juce::Justification::left);
        curX += 52.0f;

        // Swatch 4: Tonal Peaks (Only when Reduction is Unlinked)
        if (showTonalSwatch)
        {
            g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks);
            juce::Line<float> dashLine(curX, legendY + 11.0f, curX + 13.0f, legendY + 11.0f);
            float dLen[] = { 2.0f, 2.0f };
            g.drawDashedLine(dashLine, dLen, 2, 1.5f);
            curX += 17.0f;
            g.setColour(juce::Colour(0xffd8e0ec));
            g.drawText("Tonal Peaks", static_cast<int>(curX), static_cast<int>(legendY), 65, static_cast<int>(legendH), juce::Justification::left);
        }
    }
}
