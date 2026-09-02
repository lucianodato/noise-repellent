/*
noise-repellent-live -- Zero-Latency Noise Reduction JUCE Plugin

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

#include "LiveSpectralVisualizer.h"
#include "../../Shared/GUI/LookAndFeel.h"
#include <algorithm>
#include <cmath>

// ERB band center frequencies, copied from libspecbleach critical_bands.c
// (erb_bands[38]) so the display matches the processor's band layout.
static const float erbBandCenters[LiveSpectralVisualizerComponent::kNumBands] =
    {64.F,    98.F,    135.F,   176.F,   222.F,   273.F,  330.F,  393.F,
     464.F,   545.F,   634.F,   731.F,   837.F,   953.F,  1083.F, 1235.F,
     1411.F,  1613.F,  1848.F,  2062.F,  2314.F,  2612.F, 2945.F, 3321.F,
     3752.F,  4245.F,  4795.F,  5418.F,  6128.F,  6951.F, 7893.F, 8975.F,
     10220.F, 11641.F, 13257.F, 15096.F, 17196.F, 19473.F};

// ERB-rate position (Glasberg & Moore): n = 21.4 * log10(0.00437 f + 1)
static float freqToErbRate(float freqHz) {
  return 21.4f * std::log10(0.00437f * std::max(freqHz, 1.0f) + 1.0f);
}

LiveSpectralVisualizerComponent::LiveSpectralVisualizerComponent(
    NoiseRepellentLiveAudioProcessor& p)
    : processor(p) {
  smoothedInputDB.fill(-100.0f);
  smoothedOutputDB.fill(-100.0f);
  bandFloorDB.fill(-100.0f);
  setInterceptsMouseClicks(false, false);
  startTimerHz(60);
}

LiveSpectralVisualizerComponent::~LiveSpectralVisualizerComponent() {
  stopTimer();
}

void LiveSpectralVisualizerComponent::parentHierarchyChanged() {
  rebuildBandMapping();
}

void LiveSpectralVisualizerComponent::rebuildBandMapping() {
  const uint32_t sampleRate = static_cast<uint32_t>(processor.getSampleRate());
  if (sampleRate == 0 || sampleRate == sampleRateForBands) {
    return;
  }
  sampleRateForBands = sampleRate;

  const uint32_t realSpectrumSize =
      NoiseRepellentLiveAudioProcessor::kFftBins; // N/2 bins
  const float binHz =
      static_cast<float>(sampleRate) /
      static_cast<float>(NoiseRepellentLiveAudioProcessor::kFftSize);

  numValidBands = 0;
  uint32_t prevDelimiter = 1; // skip DC
  for (size_t b = 0; b < kNumBands; ++b) {
    const uint32_t delimiter =
        std::min(realSpectrumSize,
                 static_cast<uint32_t>(std::round(erbBandCenters[b] / binHz)));

    if (delimiter <= prevDelimiter) {
      continue; // band collapsed (below resolution) — skip like engine does
    }

    bandStartBin[numValidBands] = prevDelimiter;
    bandEndBin[numValidBands] = delimiter;
    ++numValidBands;
    prevDelimiter = delimiter;
  }
  if (numValidBands > 0) {
    bandEndBin[numValidBands - 1] = realSpectrumSize; // last band to Nyquist
  }

  // Mapping changed: restart learning envelope
  bandFloorDB.fill(-100.0f);
  hasLearnedFloor = false;
}

void LiveSpectralVisualizerComponent::aggregateBands(
    const std::array<float, NoiseRepellentLiveAudioProcessor::kFftBins>&
        spectrumDB,
    std::array<float, kNumBands>& bandsDB) const {
  for (uint32_t b = 0; b < numValidBands; ++b) {
    // Energy average across band bins, then back to dB
    float sum = 0.0f;
    for (uint32_t i = bandStartBin[b]; i < bandEndBin[b]; ++i) {
      sum += std::pow(10.0f, spectrumDB[i] / 10.0f);
    }
    const float count = static_cast<float>(bandEndBin[b] - bandStartBin[b]);
    bandsDB[b] = 10.0f * std::log10(std::max(sum / count, 1e-12f));
  }
}

void LiveSpectralVisualizerComponent::startLearning() {
  bandFloorDB.fill(-100.0f);
  hasLearnedFloor = false;
  learning = true;
  isSmoothedInitialized = false; // seed floor from first frame
}

void LiveSpectralVisualizerComponent::stopLearning() {
  learning = false;
}

void LiveSpectralVisualizerComponent::timerCallback() {
  if (sampleRateForBands == 0) {
    rebuildBandMapping();
  }

  const bool frameReceived = processor.getNextSpectralFrame(currentFrame);
  if (frameReceived) {
    const size_t numBins = NoiseRepellentLiveAudioProcessor::kFftBins;

    if (!isSmoothedInitialized) {
      smoothedInputDB = currentFrame.inputMagnitudeDB;
      smoothedOutputDB = currentFrame.outputMagnitudeDB;
      isSmoothedInitialized = true;
    } else {
      // Asymmetric EMA (fast attack, smooth release)
      constexpr float kAttackAlpha = 0.5f;
      constexpr float kReleaseAlpha = 0.25f;

      for (size_t i = 0; i < numBins; ++i) {
        const float targetIn = currentFrame.inputMagnitudeDB[i];
        const float alphaIn =
            (targetIn > smoothedInputDB[i]) ? kAttackAlpha : kReleaseAlpha;
        smoothedInputDB[i] += alphaIn * (targetIn - smoothedInputDB[i]);

        const float targetOut = currentFrame.outputMagnitudeDB[i];
        const float alphaOut =
            (targetOut > smoothedOutputDB[i]) ? kAttackAlpha : kReleaseAlpha;
        smoothedOutputDB[i] += alphaOut * (targetOut - smoothedOutputDB[i]);
      }
    }

    if (learning && numValidBands > 0) {
      // Learn from the raw frame (not the display-smoothed one) so the
      // envelope converges at engine speed, not display speed.
      aggregateBands(currentFrame.inputMagnitudeDB, bandInputDB);
      for (uint32_t b = 0; b < numValidBands; ++b) {
        // First frame seeds the envelope; afterwards minimum statistics
        // with a slow upward drift so a looped section converges.
        if (!hasLearnedFloor || bandInputDB[b] < bandFloorDB[b]) {
          bandFloorDB[b] = bandInputDB[b];
        } else {
          bandFloorDB[b] = std::min(bandFloorDB[b] + 0.2f, -20.0f);
        }
      }
      hasLearnedFloor = true;
    }
  }

  if (frameReceived) {
    idleTicks = 0;
    repaint();
  } else if (idleTicks < 30) {
    idleTicks++;
    repaint();
  }
}

// Map an ERB rate value to an X pixel position (linear in ERB scale)
static float erbToX(float erbRate, float width, float maxErbRate) {
  return juce::jlimit(0.0f, width, (erbRate / maxErbRate) * width);
}

void LiveSpectralVisualizerComponent::paint(juce::Graphics& g) {
  g.fillAll(juce::Colour(0xff232832));

  const float w = static_cast<float>(getWidth());
  const float h = static_cast<float>(getHeight());
  const float minDB = -100.0f;
  const float maxDB = -20.0f;
  const float maxErbRate = freqToErbRate(20000.0f);

  auto dbToY = [&](float db) {
    const float clamped = juce::jlimit(minDB, maxDB, db);
    return h * (1.0f - (clamped - minDB) / (maxDB - minDB));
  };
  auto centerErb = [&](uint32_t b) {
    const uint32_t start = bandStartBin[b];
    const uint32_t end =
        (b + 1 < numValidBands) ? bandStartBin[b + 1] : bandEndBin[b];
    const float centerBin = static_cast<float>(start + end) * 0.5f;
    const float centerHz =
        centerBin * static_cast<float>(sampleRateForBands) /
        static_cast<float>(NoiseRepellentLiveAudioProcessor::kFftSize);
    return freqToErbRate(centerHz);
  };

  // ── Grid lines ──
  g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                              juce::Font::bold));

  for (int db = -90; db <= -20; db += 10) {
    const float y = dbToY(static_cast<float>(db));
    g.setColour(juce::Colour(0xff3d4657));
    g.drawHorizontalLine(static_cast<int>(y), 0.0f, w);

    g.setColour(juce::Colour(0xffa8b3c4));
    g.drawText(juce::String(db) + " dB", 8, static_cast<int>(y) - 12, 60, 12,
               juce::Justification::left);
  }

  static const float freqLabels[] = {50,   100,  200,   500,  1000,
                                     2000, 5000, 10000, 20000};
  for (float f : freqLabels) {
    if (f >= sampleRateForBands / 2.0f) {
      continue;
    }
    const float x = erbToX(freqToErbRate(f), w, maxErbRate);
    g.setColour(juce::Colour(0xff3d4657));
    g.drawVerticalLine(static_cast<int>(x), 0.0f, h);

    g.setColour(juce::Colour(0xffa8b3c4));
    const juce::String label = f >= 1000.0f ? juce::String(f / 1000.0f, 0) + "k"
                                            : juce::String(static_cast<int>(f));
    g.drawText(label, static_cast<int>(x) + 4, static_cast<int>(h) - 14, 40, 12,
               juce::Justification::left);
  }

  if (numValidBands == 0) {
    return;
  }

  aggregateBands(smoothedInputDB, bandInputDB);
  aggregateBands(smoothedOutputDB, bandOutputDB);
  for (uint32_t b = 0; b < numValidBands; ++b) {
    bandDeltaDB[b] = bandInputDB[b] - bandOutputDB[b];
  }

  // ── Learned Noise Floor Shape (amber profile line) ──
  if (hasLearnedFloor) {
    juce::Path floorPath;
    for (uint32_t b = 0; b < numValidBands; ++b) {
      const float x = erbToX(centerErb(b), w, maxErbRate);
      const float y = dbToY(bandFloorDB[b]);
      if (b == 0) {
        floorPath.startNewSubPath(x, y);
      } else {
        floorPath.lineTo(x, y);
      }
    }
    g.setColour(
        learning ? NoiseRepellentLookAndFeel::kColorNoiseProfile.brighter(0.3f)
                 : NoiseRepellentLookAndFeel::kColorNoiseProfile);
    g.strokePath(floorPath, juce::PathStrokeType(2.0f));
  }

  // ── Input Signal Spectrum (Filled Translucent Area, ERB bands) ──
  {
    juce::Path inputAreaPath;
    float lastX = 0.0f;
    for (uint32_t b = 0; b < numValidBands; ++b) {
      const float x = erbToX(centerErb(b), w, maxErbRate);
      const float y = dbToY(bandInputDB[b]);
      if (b == 0) {
        inputAreaPath.startNewSubPath(0.0f, h);
        inputAreaPath.lineTo(0.0f, y);
        inputAreaPath.lineTo(x, y);
      } else {
        inputAreaPath.lineTo(x, y);
      }
      lastX = x;
    }
    inputAreaPath.lineTo(lastX, h);
    inputAreaPath.lineTo(0.0f, h);
    inputAreaPath.closeSubPath();

    g.setColour(NoiseRepellentLookAndFeel::kColorInputSignal.withAlpha(0.30f));
    g.fillPath(inputAreaPath);
  }

  // ── Denoised Output Signal Curve (Bright Solid Line, ERB bands) ──
  {
    juce::Path outputPath;
    for (uint32_t b = 0; b < numValidBands; ++b) {
      const float x = erbToX(centerErb(b), w, maxErbRate);
      const float y = dbToY(bandOutputDB[b]);
      if (b == 0) {
        outputPath.startNewSubPath(x, y);
      } else {
        outputPath.lineTo(x, y);
      }
    }

    g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
    g.strokePath(outputPath, juce::PathStrokeType(1.8f));
  }

  // ── Delta Curve (reduction applied per band) ──
  if (deltaVisible) {
    juce::Path deltaPath;
    for (uint32_t b = 0; b < numValidBands; ++b) {
      const float x = erbToX(centerErb(b), w, maxErbRate);
      const float y = dbToY(juce::jlimit(minDB, 0.0f, bandDeltaDB[b] + minDB));
      if (b == 0) {
        deltaPath.startNewSubPath(x, y);
      } else {
        deltaPath.lineTo(x, y);
      }
    }

    g.setColour(juce::Colour(0xffb48ead));
    g.strokePath(deltaPath, juce::PathStrokeType(2.0f));
  }

  // ── Color Legend (Top-Center Overlay) ──
  {
    constexpr float padding = 10.0f;
    const float swatch1W = 10.0f + 4.0f + 32.0f + 14.0f; // Input (60)
    const float swatch2W = 12.0f + 4.0f + 40.0f + 14.0f; // Output (70)
    const float swatch3W = 12.0f + 4.0f + 34.0f;         // Delta
    const bool showDelta = deltaVisible;
    const float legendW =
        padding + swatch1W + swatch2W + (showDelta ? swatch3W : 0.0f) + padding;
    const float legendH = 24.0f;
    const float legendX = (w - legendW) * 0.5f;
    const float legendY = 10.0f;

    if (legendX >= 10.0f && (legendX + legendW) <= w - 10.0f) {
      g.setColour(juce::Colour(0xeb252a35));
      g.fillRoundedRectangle(legendX, legendY, legendW, legendH, 4.0f);
      g.setColour(juce::Colour(0xff4c566a));
      g.drawRoundedRectangle(legendX, legendY, legendW, legendH, 4.0f, 1.0f);

      g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                                  juce::Font::bold));

      float curX = legendX + padding;

      g.setColour(
          NoiseRepellentLookAndFeel::kColorInputSignal.withAlpha(0.70f));
      g.fillRect(curX, legendY + 7.0f, 10.0f, 8.0f);
      curX += 14.0f;
      g.setColour(juce::Colour(0xffd8e0ec));
      g.drawText("Input", static_cast<int>(curX), static_cast<int>(legendY), 32,
                 static_cast<int>(legendH), juce::Justification::left);
      curX += 32.0f + 14.0f;

      g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
      g.drawLine(curX, legendY + 11.0f, curX + 12.0f, legendY + 11.0f, 2.0f);
      curX += 16.0f;
      g.setColour(juce::Colour(0xffd8e0ec));
      g.drawText("Output", static_cast<int>(curX), static_cast<int>(legendY),
                 40, static_cast<int>(legendH), juce::Justification::left);
      curX += 40.0f + 14.0f;

      if (showDelta) {
        g.setColour(juce::Colour(0xffb48ead));
        g.drawLine(curX, legendY + 11.0f, curX + 12.0f, legendY + 11.0f, 2.0f);
        curX += 16.0f;
        g.setColour(juce::Colour(0xffd8e0ec));
        g.drawText("Delta", static_cast<int>(curX), static_cast<int>(legendY),
                   34, static_cast<int>(legendH), juce::Justification::left);
      }
    }
  }
}
