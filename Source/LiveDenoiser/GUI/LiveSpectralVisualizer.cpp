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

// Bark-rate position (Zwicker): b = 13*atan(0.76f/1000) + 3.5*atan((f/7500)^2)
static float freqToBarkRate(float freqHz) {
  const float f = std::max(freqHz, 1.0f);
  return 13.0f * std::atan(0.76f * f / 1000.0f) +
         3.5f * std::atan((f / 7500.0f) * (f / 7500.0f));
}

static float linearToDB(float linear) {
  return 20.0f * std::log10(std::max(linear, 1e-7f));
}

LiveSpectralVisualizerComponent::LiveSpectralVisualizerComponent(
    NoiseRepellentLiveAudioProcessor& p)
    : processor(p) {
  smoothedInputDB.fill(-100.0f);
  smoothedOutputDB.fill(-100.0f);
  smoothedThresholdDB.fill(-100.0f);
  setInterceptsMouseClicks(false, false);
  startTimerHz(60);
}

LiveSpectralVisualizerComponent::~LiveSpectralVisualizerComponent() {
  stopTimer();
}

void LiveSpectralVisualizerComponent::timerCallback() {
  // Refresh the engine scale (cheap copy; defines the axis span).
  if (processor.getLiveBandEdges(bandLoHz.data(), bandHiHz.data())) {
    numActiveBands = 0;
    for (size_t b = 0; b < kNumBands; ++b) {
      if (bandHiHz[b] <= 0.0f) {
        break;
      }
      ++numActiveBands;
    }
    if (numActiveBands > 0) {
      axisBarkLo = freqToBarkRate(bandLoHz[0]);
      axisBarkHi = freqToBarkRate(bandHiHz[numActiveBands - 1]);
      if (!(axisBarkHi > axisBarkLo)) {
        axisBarkHi = axisBarkLo + 1.0f;
      }
    }
  }

  // Drain to the latest frame so the display never lags behind.
  NoiseRepellentLiveAudioProcessor::BandFrame frame;
  bool frameReceived = false;
  while (processor.getNextBandFrame(frame)) {
    frameReceived = true;
  }
  if (!frameReceived) {
    if (idleTicks < 30) {
      idleTicks++;
      repaint();
    }
    return;
  }
  idleTicks = 0;

  // Asymmetric EMA (fast attack, smooth release) on the display values.
  constexpr float kAttackAlpha = 0.5f;
  constexpr float kReleaseAlpha = 0.25f;

  for (size_t b = 0; b < kNumBands; ++b) {
    const float targets[3] = {linearToDB(frame.inputLevels[b]),
                              linearToDB(frame.outputLevels[b]),
                              linearToDB(frame.thresholdLevels[b])};
    float* smoothed[3] = {&smoothedInputDB[b], &smoothedOutputDB[b],
                          &smoothedThresholdDB[b]};
    for (int c = 0; c < 3; ++c) {
      if (!isSmoothedInitialized) {
        *smoothed[c] = targets[c];
      } else {
        const float alpha =
            (targets[c] > *smoothed[c]) ? kAttackAlpha : kReleaseAlpha;
        *smoothed[c] += alpha * (targets[c] - *smoothed[c]);
      }
    }
  }
  isSmoothedInitialized = true;
  repaint();
}

void LiveSpectralVisualizerComponent::paint(juce::Graphics& g) {
  g.fillAll(juce::Colour(0xff232832));

  const float w = static_cast<float>(getWidth());
  const float h = static_cast<float>(getHeight());
  const float minDB = -100.0f;
  const float maxDB = 0.0f;

  // Plot area: dB scale on the right (RX style), Hz labels along the
  // bottom. The X axis spans exactly the engine's active Bark range so
  // the curves fill the display at any sample rate.
  constexpr float kRightMargin = 46.0f;
  constexpr float kBottomMargin = 16.0f;
  constexpr float kTopMargin = 14.0f; // keeps the 0 dB label off the border
  const float plotW = std::max(w - kRightMargin - 4.0f, 10.0f);
  const float plotH = std::max(h - kTopMargin - kBottomMargin - 4.0f, 10.0f);
  const float barkSpan = std::max(axisBarkHi - axisBarkLo, 1e-6f);

  auto dbToY = [&](float db) {
    const float clamped = juce::jlimit(minDB, maxDB, db);
    return kTopMargin + plotH * (1.0f - (clamped - minDB) / (maxDB - minDB));
  };
  auto barkToX = [&](float barkRate) {
    return 2.0f + juce::jlimit(0.0f, plotW,
                               ((barkRate - axisBarkLo) / barkSpan) * plotW);
  };
  auto centerBark = [&](size_t b) {
    const float centerHz =
        std::sqrt(std::max(bandLoHz[b], 1.0f) * std::max(bandHiHz[b], 2.0f));
    return freqToBarkRate(centerHz);
  };

  // ── Grid lines ──
  g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                              juce::Font::bold));

  for (int db = -90; db <= 0; db += 10) {
    const float y = dbToY(static_cast<float>(db));
    g.setColour(juce::Colour(0xff3d4657));
    g.drawHorizontalLine(static_cast<int>(y), 0.0f, w);

    g.setColour(juce::Colour(0xffa8b3c4));
    g.drawText(juce::String(db) + " dB", static_cast<int>(w) - 44,
               static_cast<int>(y) - 6, 42, 12,
               juce::Justification::centredRight);
  }

  static const float freqLabels[] = {100,  200,  500,   1000,
                                     2000, 5000, 10000, 20000};
  for (float f : freqLabels) {
    const float bark = freqToBarkRate(f);
    if (bark < axisBarkLo || bark > axisBarkHi) {
      continue;
    }
    const float x = barkToX(bark);
    g.setColour(juce::Colour(0xff3d4657));
    g.drawVerticalLine(static_cast<int>(x), 0.0f, h - kBottomMargin + 4.0f);

    g.setColour(juce::Colour(0xffa8b3c4));
    const juce::String label = f >= 1000.0f ? juce::String(f / 1000.0f, 0) + "k"
                                            : juce::String(static_cast<int>(f));
    g.drawText(label, static_cast<int>(x) - 20,
               static_cast<int>(h) - static_cast<int>(kBottomMargin), 40, 12,
               juce::Justification::centred);
  }

  if (numActiveBands == 0) {
    return;
  }

  // ── Gate Threshold Shape (amber profile line, always visible) ──
  {
    juce::Path thresholdPath;
    for (size_t b = 0; b < numActiveBands; ++b) {
      const float x = barkToX(centerBark(b));
      const float y = dbToY(smoothedThresholdDB[b]);
      if (b == 0) {
        thresholdPath.startNewSubPath(x, y);
      } else {
        thresholdPath.lineTo(x, y);
      }
    }
    g.setColour(NoiseRepellentLookAndFeel::kColorNoiseProfile);
    g.strokePath(thresholdPath, juce::PathStrokeType(2.0f));
  }

  // ── Input Signal Energy (Filled Translucent Area, Bark bands) ──
  {
    juce::Path inputAreaPath;
    float lastX = 0.0f;
    for (size_t b = 0; b < numActiveBands; ++b) {
      const float x = barkToX(centerBark(b));
      const float y = dbToY(smoothedInputDB[b]);
      if (b == 0) {
        inputAreaPath.startNewSubPath(x, dbToY(minDB));
        inputAreaPath.lineTo(x, y);
      } else {
        inputAreaPath.lineTo(x, y);
      }
      lastX = x;
    }
    inputAreaPath.lineTo(lastX, dbToY(minDB));
    inputAreaPath.closeSubPath();

    g.setColour(NoiseRepellentLookAndFeel::kColorInputSignal.withAlpha(0.30f));
    g.fillPath(inputAreaPath);
  }

  // ── Denoised Output Energy Curve (Bright Solid Line, Bark bands) ──
  {
    juce::Path outputPath;
    for (size_t b = 0; b < numActiveBands; ++b) {
      const float x = barkToX(centerBark(b));
      const float y = dbToY(smoothedOutputDB[b]);
      if (b == 0) {
        outputPath.startNewSubPath(x, y);
      } else {
        outputPath.lineTo(x, y);
      }
    }

    g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
    g.strokePath(outputPath, juce::PathStrokeType(1.8f));
  }

  // ── Color Legend (Top-Center Overlay) ──
  {
    constexpr float padding = 10.0f;
    const float inputW = 10.0f + 4.0f + 32.0f + 14.0f;     // Input (60)
    const float outputW = 12.0f + 4.0f + 40.0f + 14.0f;    // Output (70)
    const float thresholdW = 12.0f + 4.0f + 56.0f + 14.0f; // Threshold (96)
    const float legendW = padding + inputW + outputW + thresholdW + padding;
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

      g.setColour(NoiseRepellentLookAndFeel::kColorNoiseProfile);
      g.drawLine(curX, legendY + 11.0f, curX + 12.0f, legendY + 11.0f, 2.0f);
      curX += 16.0f;
      g.setColour(juce::Colour(0xffd8e0ec));
      g.drawText("Threshold", static_cast<int>(curX), static_cast<int>(legendY),
                 56, static_cast<int>(legendH), juce::Justification::left);
    }
  }
}