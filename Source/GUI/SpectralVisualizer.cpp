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
#include <algorithm>
#include <cmath>

namespace {

// Tooltip strings (single source; getTooltip() only selects among these)
const juce::String kTpBadgeTip =
    "Transient Protection (TP): Click to toggle.\n"
    "Preserves sharp attack transients and plucked notes.";

} // namespace

SpectralVisualizerComponent::SpectralVisualizerComponent(
    NoiseRepellentAudioProcessor& p)
    : processor(p) {
  smoothedInputDB.fill(SpectralVisualizerComponent::kAxisMinDB);
  smoothedOutputDB.fill(SpectralVisualizerComponent::kAxisMinDB);
  startTimerHz(60);
}

SpectralVisualizerComponent::~SpectralVisualizerComponent() {
  stopTimer();
}

void SpectralVisualizerComponent::timerCallback() {
  bool frameReceived = processor.getNextSpectralFrame(currentFrame);
  if (frameReceived) {
    const size_t numBins = NoiseRepellentAudioProcessor::kFftBins;

    if (!isSmoothedInitialized) {
      smoothedInputDB = currentFrame.inputMagnitudeDB;
      smoothedOutputDB = currentFrame.outputMagnitudeDB;
      isSmoothedInitialized = true;
    } else {
      // Asymmetric Exponential Moving Average (Fast Attack, Smooth Fluid Decay)
      constexpr float kAttackAlpha = 0.35f; // Fast response to audio transients
      constexpr float kReleaseAlpha =
          0.10f; // Smooth, liquid decay on falling signals

      for (size_t i = 0; i < numBins; ++i) {
        float targetIn = currentFrame.inputMagnitudeDB[i];
        float alphaIn =
            (targetIn > smoothedInputDB[i]) ? kAttackAlpha : kReleaseAlpha;
        smoothedInputDB[i] += alphaIn * (targetIn - smoothedInputDB[i]);

        float targetOut = currentFrame.outputMagnitudeDB[i];
        float alphaOut =
            (targetOut > smoothedOutputDB[i]) ? kAttackAlpha : kReleaseAlpha;
        smoothedOutputDB[i] += alphaOut * (targetOut - smoothedOutputDB[i]);
      }
    }
  }

  // Synchronize transient detection active state and LED activity
  auto* transientParam =
      processor.getAPVTS().getRawParameterValue("transient_protection_enable");
  transientProtectionActive =
      (transientParam != nullptr && transientParam->load() > 0.5f);

  if (!transientProtectionActive) {
    ledBrightness = 0.0f;
    transientHoldTicks = 0;
  } else {
    // Synchronize transient detection LED strictly with the displayed
    // visualizer FFT frame
    float targetIntensity = currentFrame.transientIntensity;

    // Trigger LED strictly when high-confidence transient is present in the
    // displayed frame
    if (targetIntensity >= 0.35f) {
      transientHoldTicks = 4; // ~65ms crisp hold at 60Hz
      ledBrightness = 1.0f;   // Full crisp light
    } else if (transientHoldTicks > 0) {
      transientHoldTicks--;
    } else {
      ledBrightness = 0.0f; // Immediate clean turn-off
    }
  }

  if (frameReceived || transientHoldTicks > 0 || ledBrightness > 0.0f) {
    idleTicks = 0;
    repaint();
  } else if (idleTicks < 30) {
    idleTicks++;
    repaint();
  }
}

void SpectralVisualizerComponent::resized() {
}

juce::Rectangle<float> SpectralVisualizerComponent::getTpBadgeBounds() const {
  const float w = static_cast<float>(getWidth());
  return {w - kBadgeW - kBadgeMargin, kBadgeTop, kBadgeW, kBadgeH};
}

float SpectralVisualizerComponent::biasToY(float biasDB, float h) {
  return h * 0.5f - (biasDB / kCurveMaxBiasDB) * (h * kCurveHeightFrac);
}

float SpectralVisualizerComponent::yToBias(float y, float h) {
  return ((h * 0.5f - y) / (h * kCurveHeightFrac)) * kCurveMaxBiasDB;
}

// Map a frequency (Hz) to an X pixel position using logarithmic scale
static float freqToX(float freqHz, float width,
                     float minFreq = SpectralVisualizerComponent::kAxisMinFreq,
                     float maxFreq = SpectralVisualizerComponent::kAxisMaxFreq) {
  if (freqHz <= minFreq)
    return 0.0f;
  if (freqHz >= maxFreq)
    return width;
  return (std::log10(freqHz / minFreq) / std::log10(maxFreq / minFreq)) * width;
}

// Map a dB value to a Y pixel position
static float dbToY(float db, float height,
                   float minDB = SpectralVisualizerComponent::kAxisMinDB,
                   float maxDB = SpectralVisualizerComponent::kAxisMaxDB) {
  float clamped = juce::jlimit(minDB, maxDB, db);
  float normalized =
      (clamped - minDB) / (maxDB - minDB); // 0 at bottom, 1 at top
  return height * (1.0f - normalized);
}

// Apply variable fractional-octave spatial frequency smoothing across FFT bins
// (UI display only). Low frequencies (bass/hum) get minimal smoothing
// (preserving sharp peaks), while high frequencies get progressively wider
// smoothing windows matching log-scale visuals.
static void smoothBinsLogarithmicFrequencyDomain(const float* src, float* dst,
                                                 size_t numBins) {
  if (numBins == 0)
    return;

  for (size_t i = 0; i < numBins; ++i) {
    // Half-window radius grows proportionally with bin index i (fractional
    // octave) i = 5 (~50 Hz)   -> radius = 0 (exact bin value) i = 20 (~200 Hz)
    // -> radius = 1 (3-bin window) i = 100 (~1 kHz) -> radius = 3 (7-bin
    // window) i = 500 (~5 kHz) -> radius = 12 (25-bin window) i = 1000 (~10
    // kHz)-> radius = 20 (41-bin window)
    int radius = std::clamp(
        static_cast<int>(std::round(static_cast<float>(i) * 0.025f)), 0, 20);

    if (radius == 0) {
      dst[i] = src[i];
      continue;
    }

    int minIdx = std::max(0, static_cast<int>(i) - radius);
    int maxIdx =
        std::min(static_cast<int>(numBins - 1), static_cast<int>(i) + radius);

    float sum = 0.0f;
    float totalWeight = 0.0f;

    for (int k = minIdx; k <= maxIdx; ++k) {
      // Triangular window weight: highest weight at center bin i, tapering down
      // to edges
      float dist = std::abs(static_cast<float>(k - static_cast<int>(i)));
      float weight = 1.0f - (dist / static_cast<float>(radius + 1));

      sum += src[k] * weight;
      totalWeight += weight;
    }

    dst[i] = (totalWeight > 0.0f) ? (sum / totalWeight) : src[i];
  }
}

void SpectralVisualizerComponent::paint(juce::Graphics& g) {
  g.fillAll(juce::Colour(0xff232832));

  const float w = static_cast<float>(getWidth());
  const float h = static_cast<float>(getHeight());
  const float minDB = SpectralVisualizerComponent::kAxisMinDB;
  const float maxDB = SpectralVisualizerComponent::kAxisMaxDB;
  const float minFreq = SpectralVisualizerComponent::kAxisMinFreq;
  const float maxFreq = SpectralVisualizerComponent::kAxisMaxFreq;

  const double sampleRate = processor.getSampleRate();
  const size_t numBins = NoiseRepellentAudioProcessor::kFftBins;
  const float binWidth =
      static_cast<float>(sampleRate) /
      static_cast<float>(NoiseRepellentAudioProcessor::kFftSize);

  // Frequency-domain spatially smoothed display buffers (UI display smoothing
  // only). The noise profile is deliberately NOT smoothed: it arrives at the
  // engine's native bin resolution (nearest-bin mapped) and smoothing would
  // hide the coarser/finer steps that distinguish frame sizes.
  std::vector<float> freqSmoothedInput(numBins);
  std::vector<float> freqSmoothedOutput(numBins);

  smoothBinsLogarithmicFrequencyDomain(smoothedInputDB.data(),
                                       freqSmoothedInput.data(), numBins);
  smoothBinsLogarithmicFrequencyDomain(smoothedOutputDB.data(),
                                       freqSmoothedOutput.data(), numBins);

  // ── Grid lines ──
  g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                              juce::Font::bold));

  // dB Y-grid (-90 dB to -20 dB in 10 dB steps)
  for (int db = -90; db <= -20; db += 10) {
    float y = dbToY(static_cast<float>(db), h, minDB, maxDB);
    g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorGridLine));
    g.drawHorizontalLine(static_cast<int>(y), 0.0f, w);

    g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorGridLabel));
    juce::String label = juce::String(db) + " dB";
    g.drawText(label, 8, static_cast<int>(y) - 12, 60, 12,
               juce::Justification::left);
  }

  // Frequency X-grid (logarithmic)
  static const float freqs[] = {50,   100,  200,   500,  1000,
                                2000, 5000, 10000, 20000};
  for (float f : freqs) {
    float x = freqToX(f, w, minFreq, maxFreq);
    g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorGridLine));
    g.drawVerticalLine(static_cast<int>(x), 0.0f, h);

    g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorGridLabel));
    juce::String label = f >= 1000.0f ? juce::String(f / 1000.0f, 0) + "k"
                                      : juce::String(static_cast<int>(f));
    g.drawText(label, static_cast<int>(x) + 4, static_cast<int>(h) - 14, 40, 12,
               juce::Justification::left);
  }

  // ── Helper: build a path by mapping FFT bins to log-frequency X positions ──
  auto buildLogFreqPath = [&](const float* dbData, size_t bins) -> juce::Path {
    juce::Path path;
    bool started = false;

    for (size_t i = 1; i < bins; ++i) // skip bin 0 (DC)
    {
      float freq = static_cast<float>(i) * binWidth;
      if (freq < minFreq || freq > maxFreq)
        continue;

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

    for (size_t i = 1; i < numBins; ++i) {
      float freq = static_cast<float>(i) * binWidth;
      if (freq < minFreq || freq > maxFreq)
        continue;

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

      g.setColour(
          NoiseRepellentLookAndFeel::kColorInputSignal.withAlpha(0.30f));
      g.fillPath(inputAreaPath);
    }
  }

  // 2. Denoised Output Signal Curve (Bright Solid Cyan Line)
  {
    juce::Path outputPath =
        buildLogFreqPath(freqSmoothedOutput.data(), numBins);
    g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
    g.strokePath(outputPath, juce::PathStrokeType(1.8f));
  }

  // 3. Noise Floor Profile Curve (Solid Warm Amber Line, native engine bins)
  if (currentFrame.hasNoiseProfile) {
    juce::Path noisePath =
        buildLogFreqPath(currentFrame.noiseFloorDB.data(), numBins);
    g.setColour(NoiseRepellentLookAndFeel::kColorNoiseProfile);
    g.strokePath(noisePath, juce::PathStrokeType(2.0f));

    // When Threshold Offset is unlinked, redraw the profile segments around
    // detected tonal peaks using the tonal colour so the tonal threshold
    // offset is visually distinguished from the broadband noise profile.
    if (isAdvancedVisible && !currentFrame.isOffsetLinked &&
        !currentFrame.tonalPeaksHz.empty()) {
      constexpr float kTonalRegionRatio = 0.05f;

      std::vector<bool> tonalBin(numBins, false);
      for (size_t i = 1; i < numBins; ++i) {
        float freq = static_cast<float>(i) * binWidth;
        for (float peakHz : currentFrame.tonalPeaksHz) {
          if (freq >= peakHz * (1.0f - kTonalRegionRatio) &&
              freq <= peakHz * (1.0f + kTonalRegionRatio)) {
            tonalBin[i] = true;
            break;
          }
        }
      }

      size_t i = 1;
      while (i < numBins) {
        if (!tonalBin[i]) {
          ++i;
          continue;
        }

        size_t runEnd = i;
        while (runEnd + 1 < numBins && tonalBin[runEnd + 1])
          ++runEnd;

        juce::Path segment;
        bool started = false;
        for (size_t b = i; b <= runEnd + 1 && b < numBins; ++b) {
          float freq = static_cast<float>(b) * binWidth;
          if (freq < minFreq || freq > maxFreq)
            continue;
          float x = freqToX(freq, w, minFreq, maxFreq);
          float y = dbToY(currentFrame.noiseFloorDB[b], h, minDB, maxDB);
          if (!started) {
            segment.startNewSubPath(x, y);
            started = true;
          } else {
            segment.lineTo(x, y);
          }
        }
        if (started) {
          g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks);
          g.strokePath(segment, juce::PathStrokeType(2.6f));
        }

        i = runEnd + 1;
      }
    }
  }

  // 3b. Reduction Curve Bias Overlay & Nodes (Only shown in Advanced mode)
  if (isAdvancedVisible && currentFrame.reductionCurveEnabled) {
    const auto& nodes = processor.getCurveNodes();
    if (nodes.size() >= 2) {
      auto nodeToPoint = [&](const NoiseRepellentAudioProcessor::CurveNode& n)
          -> juce::Point<float> {
        float nx = std::clamp(n.normX, 0.0f, 1.0f);
        float ny = SpectralVisualizerComponent::biasToY(n.biasDB, h);
        return {nx * w, ny};
      };

      // Draw zero-bias centerline (dashed soft green)
      g.setColour(
          NoiseRepellentLookAndFeel::kColorReductionCurve.withAlpha(0.35f));
      juce::Line<float> zeroLine(0.0f, h * 0.5f, w, h * 0.5f);
      float dashLen[] = {4.0f, 4.0f};
      g.drawDashedLine(zeroLine, dashLen, 2, 1.0f);

      // Draw smooth cubic spline curve path using JUCE Path cubicTo
      juce::Path curvePath;
      juce::Point<float> pStart = nodeToPoint(nodes.front());
      curvePath.startNewSubPath(pStart);

      size_t numNodes = nodes.size();
      if (numNodes == 2) {
        curvePath.lineTo(nodeToPoint(nodes.back()));
      } else {
        std::vector<juce::Point<float>> pts(numNodes);
        for (size_t i = 0; i < numNodes; ++i)
          pts[i] = nodeToPoint(nodes[i]);

        std::vector<float> dY(numNodes - 1);
        std::vector<float> dX(numNodes - 1);
        for (size_t i = 0; i < numNodes - 1; ++i) {
          dX[i] = pts[i + 1].x - pts[i].x;
          dY[i] = pts[i + 1].y - pts[i].y;
        }

        std::vector<float> m(numNodes, 0.0f);
        m.front() = (dX.front() > 1e-5f) ? dY.front() / dX.front() : 0.0f;
        m.back() = (dX.back() > 1e-5f) ? dY.back() / dX.back() : 0.0f;

        for (size_t i = 1; i < numNodes - 1; ++i) {
          float secant1 = (dX[i - 1] > 1e-5f) ? dY[i - 1] / dX[i - 1] : 0.0f;
          float secant2 = (dX[i] > 1e-5f) ? dY[i] / dX[i] : 0.0f;
          if (secant1 * secant2 <= 0.0f)
            m[i] = 0.0f;
          else
            m[i] = 0.5f * (secant1 + secant2);
        }

        for (size_t i = 0; i < numNodes - 1; ++i) {
          float dx = dX[i];
          juce::Point<float> c1(pts[i].x + (1.0f / 3.0f) * dx,
                                pts[i].y + (1.0f / 3.0f) * dx * m[i]);
          juce::Point<float> c2(pts[i + 1].x - (1.0f / 3.0f) * dx,
                                pts[i + 1].y - (1.0f / 3.0f) * dx * m[i + 1]);
          curvePath.cubicTo(c1, c2, pts[i + 1]);
        }
      }

      g.setColour(NoiseRepellentLookAndFeel::kColorReductionCurve);
      g.strokePath(curvePath, juce::PathStrokeType(2.2f));

      // Draw nodes
      for (size_t i = 0; i < nodes.size(); ++i) {
        juce::Point<float> pt = nodeToPoint(nodes[i]);
        g.setColour(NoiseRepellentLookAndFeel::kColorReductionCurve);
        g.fillEllipse(pt.x - 5.0f, pt.y - 5.0f, 10.0f, 10.0f);
        g.setColour(juce::Colours::white);
        g.drawEllipse(pt.x - 5.0f, pt.y - 5.0f, 10.0f, 10.0f, 1.5f);
      }

      // Draw dB Reference Y-Scale for Reduction Curve on Right Margin
      const float scaleX = w - 42.0f;
      const float textX = w - 38.0f;
      static const int biasLevels[] = {+24, +12, 0, -12, -24}; // +/-kCurveMaxBiasDB

      for (int biasVal : biasLevels) {
        float ny =
            SpectralVisualizerComponent::biasToY(static_cast<float>(biasVal), h);
        g.setColour(
            NoiseRepellentLookAndFeel::kColorReductionCurve.withAlpha(0.35f));
        g.drawHorizontalLine(static_cast<int>(ny), scaleX, w);

        g.setColour(NoiseRepellentLookAndFeel::kColorReductionCurve);
        g.setFont(juce::FontOptions(10.0f, juce::Font::bold));
        int reductionVal = biasVal;
        juce::String labelStr =
            (reductionVal > 0 ? "+" : "") + juce::String(reductionVal) + "dB";
        g.drawText(labelStr, static_cast<int>(textX), static_cast<int>(ny) - 7,
                   36, 14, juce::Justification::left);
      }
    }
  }

  // 4. Tonal Peak Markers (Dashed vertical lines with staggered multi-tier
  // frequency tags; shown ONLY in Advanced mode when Reduction or Threshold
  // Offset is UNLINKED)
  if (isAdvancedVisible &&
      (!currentFrame.isLinked || !currentFrame.isOffsetLinked) &&
      currentFrame.hasNoiseProfile && !currentFrame.tonalPeaksHz.empty()) {
    // 3 Y-tiers starting at Y = 42.0f to guarantee zero collision with top HUD
    // overlay bar
    constexpr float kStartTagY = 42.0f;
    float lastTagRight[3] = {-100.0f, -100.0f, -100.0f};

    for (float peakHz : currentFrame.tonalPeaksHz) {
      float x = freqToX(peakHz, w, minFreq, maxFreq);
      if (x <= 4.0f || x >= w - 4.0f)
        continue;

      // Dashed vertical line across the display
      juce::Line<float> line(x, 0.0f, x, h);
      float dashLengths[] = {3.0f, 3.0f};
      g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks.withAlpha(0.75f));
      g.drawDashedLine(line, dashLengths, 2, 1.5f);

      // Frequency Tag Badge
      juce::String tag = peakHz >= 1000.0f
                             ? juce::String(peakHz / 1000.0f, 1) + "k"
                             : juce::String(static_cast<int>(peakHz)) + "Hz";
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

      if (chosenTier != -1) {
        float tagY = kStartTagY + static_cast<float>(chosenTier) * 16.0f;

        // Draw vertical leader connector line from top to badge
        g.setColour(
            NoiseRepellentLookAndFeel::kColorTonalPeaks.withAlpha(0.5f));
        g.drawVerticalLine(static_cast<int>(x), 0.0f, tagY);

        g.setColour(
            NoiseRepellentLookAndFeel::kColorTonalPeaks.withAlpha(0.9f));
        g.fillRoundedRectangle(tagX, tagY, tagW, tagH, 3.0f);

        g.setColour(juce::Colours::white);
        g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                                    juce::Font::bold));
        g.drawText(tag, static_cast<int>(tagX), static_cast<int>(tagY),
                   static_cast<int>(tagW), static_cast<int>(tagH),
                   juce::Justification::centred);

        lastTagRight[chosenTier] = tagX + tagW;
      }
    }
  }

  // 5. HUD Color Legend (Top-Right Overlay, positioned to the left of Advanced
  // Controls button)
  {
    const bool showTonalSwatch =
        isAdvancedVisible &&
        (!currentFrame.isLinked || !currentFrame.isOffsetLinked);
    const bool showCurveSwatch =
        isAdvancedVisible && currentFrame.reductionCurveEnabled;
    const float padding = 10.0f;
    const float swatch1W = 10.0f + 4.0f + 32.0f + 14.0f; // Input (60)
    const float swatch2W = 12.0f + 4.0f + 38.0f + 14.0f; // Profile (68)
    const float swatch3W =
        12.0f + 4.0f + 40.0f +
        ((showTonalSwatch || showCurveSwatch) ? 14.0f : 0.0f); // Output
    const float swatch4W =
        showTonalSwatch
            ? (12.0f + 4.0f + 64.0f + (showCurveSwatch ? 14.0f : 0.0f))
            : 0.0f; // Tonal Peaks
    const float swatch5W =
        showCurveSwatch ? (12.0f + 4.0f + 38.0f) : 0.0f; // Curve (54)

    const float legendW = padding + swatch1W + swatch2W + swatch3W + swatch4W +
                          swatch5W + padding;
    const float legendH = 24.0f;
    const float legendX = (w - legendW) * 0.5f;
    const float legendY = 10.0f;

    if (legendX >= 10.0f && (legendX + legendW) <= w - 10.0f) {
      // Semi-transparent dark glass background
      g.setColour(juce::Colour(0xeb252a35));
      g.fillRoundedRectangle(legendX, legendY, legendW, legendH, 4.0f);
      g.setColour(juce::Colour(0xff4c566a));
      g.drawRoundedRectangle(legendX, legendY, legendW, legendH, 4.0f, 1.0f);

      g.setFont(juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel,
                                  juce::Font::bold));

      float curX = legendX + padding;

      // Swatch 1: Live Input (Filled Area)
      g.setColour(
          NoiseRepellentLookAndFeel::kColorInputSignal.withAlpha(0.70f));
      g.fillRect(curX, legendY + 7.0f, 10.0f, 8.0f);
      curX += 14.0f;
      g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorLegendText));
      g.drawText("Input", static_cast<int>(curX), static_cast<int>(legendY), 32,
                 static_cast<int>(legendH), juce::Justification::left);
      curX += 32.0f + 14.0f;

      // Swatch 2: Noise Profile (Solid Amber Line)
      g.setColour(NoiseRepellentLookAndFeel::kColorNoiseProfile);
      g.drawLine(curX, legendY + 11.0f, curX + 12.0f, legendY + 11.0f, 2.0f);
      curX += 16.0f;
      g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorLegendText));
      g.drawText("Profile", static_cast<int>(curX), static_cast<int>(legendY),
                 38, static_cast<int>(legendH), juce::Justification::left);
      curX += 38.0f + 14.0f;

      // Swatch 3: Output Area (Solid Cyan Line)
      g.setColour(NoiseRepellentLookAndFeel::kColorDenoising);
      g.drawLine(curX, legendY + 11.0f, curX + 12.0f, legendY + 11.0f, 2.0f);
      curX += 16.0f;
      g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorLegendText));
      g.drawText("Output", static_cast<int>(curX), static_cast<int>(legendY),
                 40, static_cast<int>(legendH), juce::Justification::left);
      curX += 40.0f;

      // Swatch 4: Tonal Peaks (Only when Reduction is Unlinked)
      if (showTonalSwatch) {
        curX += 14.0f;
        g.setColour(NoiseRepellentLookAndFeel::kColorTonalPeaks);
        juce::Line<float> dashLine(curX, legendY + 11.0f, curX + 12.0f,
                                   legendY + 11.0f);
        float dLen[] = {2.0f, 2.0f};
        g.drawDashedLine(dashLine, dLen, 2, 1.5f);
        curX += 16.0f;
        g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorLegendText));
        g.drawText("Tonal Peaks", static_cast<int>(curX),
                   static_cast<int>(legendY), 64, static_cast<int>(legendH),
                   juce::Justification::left);
        curX += 64.0f;
      }

      // Swatch 5: Reduction Curve (Only when enabled)
      if (showCurveSwatch) {
        curX += 14.0f;
        g.setColour(NoiseRepellentLookAndFeel::kColorReductionCurve);
        g.drawLine(curX, legendY + 11.0f, curX + 12.0f, legendY + 11.0f, 2.0f);
        curX += 16.0f;
        g.setColour(juce::Colour(NoiseRepellentLookAndFeel::kColorLegendText));
        g.drawText("Curve", static_cast<int>(curX), static_cast<int>(legendY),
                   38, static_cast<int>(legendH), juce::Justification::left);
      }
    }
  }

  // 6. Transient Protection (TP) LED Indicator & Button (Top-Right Corner:
  // click to toggle ON/OFF, only visible in Advanced mode)
  if (isAdvancedVisible) {
    const juce::Rectangle<float> badge = getTpBadgeBounds();
    const float badgeX = badge.getX();
    const float badgeY = badge.getY();
    const float badgeW = badge.getWidth();
    const float badgeH = badge.getHeight();

    if (badgeX > 10.0f) {
      // Container background (dimmer when disabled)
      g.setColour(transientProtectionActive ? juce::Colour(0xeb1a1e26)
                                            : juce::Colour(0xeb14161c));
      g.fillRoundedRectangle(badgeX, badgeY, badgeW, badgeH, 4.0f);
      g.setColour(transientProtectionActive ? juce::Colour(0xff2d3542)
                                            : juce::Colour(0xff20242c));
      g.drawRoundedRectangle(badgeX, badgeY, badgeW, badgeH, 4.0f, 1.0f);

      // LED circle on the left of badge
      const float ledRadius = 3.5f;
      const float ledCenterX = badgeX + 11.0f;
      const float ledCenterY = badgeY + (badgeH * 0.5f);

      if (transientProtectionActive) {
        // Map ledBrightness to LED glow and color
        // Full bright for broadband transients (1.0), soft glow for localized
        // clicks (0.2 - 0.5) Completely off (dark) when ledBrightness == 0.0f
        juce::Colour activeLedColour =
            juce::Colour(0xffff9933); // Amber-Orange transient light

        if (ledBrightness > 0.05f) {
          // Outer glow
          g.setColour(activeLedColour.withAlpha(0.40f * ledBrightness));
          g.fillEllipse(ledCenterX - ledRadius - 3.0f,
                        ledCenterY - ledRadius - 3.0f,
                        (ledRadius + 3.0f) * 2.0f, (ledRadius + 3.0f) * 2.0f);

          juce::Colour ledFill =
              activeLedColour.withAlpha(0.20f + 0.80f * ledBrightness);
          g.setColour(ledFill);
          g.fillEllipse(ledCenterX - ledRadius, ledCenterY - ledRadius,
                        ledRadius * 2.0f, ledRadius * 2.0f);

          g.setColour(activeLedColour.withAlpha(0.40f + 0.60f * ledBrightness));
          g.drawEllipse(ledCenterX - ledRadius, ledCenterY - ledRadius,
                        ledRadius * 2.0f, ledRadius * 2.0f, 1.0f);
        } else {
          // Off state while TP is active but no transient is present
          g.setColour(juce::Colour(0xff2d3542));
          g.fillEllipse(ledCenterX - ledRadius, ledCenterY - ledRadius,
                        ledRadius * 2.0f, ledRadius * 2.0f);
          g.setColour(juce::Colour(0xff3d4756));
          g.drawEllipse(ledCenterX - ledRadius, ledCenterY - ledRadius,
                        ledRadius * 2.0f, ledRadius * 2.0f, 1.0f);
        }

        // Text label: TP ON
        g.setFont(
            juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel - 0.5f,
                              juce::Font::bold));
        g.setColour(
            juce::Colour(NoiseRepellentLookAndFeel::kColorLegendText).withAlpha(0.65f + 0.35f * ledBrightness));
        g.drawText("TP ON", static_cast<int>(badgeX + 19.0f),
                   static_cast<int>(badgeY), static_cast<int>(badgeW - 22.0f),
                   static_cast<int>(badgeH), juce::Justification::centred);
      } else {
        // Disabled state: dark off LED
        g.setColour(juce::Colour(0xff333a46));
        g.fillEllipse(ledCenterX - ledRadius, ledCenterY - ledRadius,
                      ledRadius * 2.0f, ledRadius * 2.0f);
        g.setColour(juce::Colour(0xff454f5e));
        g.drawEllipse(ledCenterX - ledRadius, ledCenterY - ledRadius,
                      ledRadius * 2.0f, ledRadius * 2.0f, 1.0f);

        g.setFont(
            juce::FontOptions(NoiseRepellentLookAndFeel::kFontSizeLabel - 0.5f,
                              juce::Font::bold));
        g.setColour(juce::Colour(0xff667080));
        g.drawText("TP OFF", static_cast<int>(badgeX + 19.0f),
                   static_cast<int>(badgeY), static_cast<int>(badgeW - 22.0f),
                   static_cast<int>(badgeH), juce::Justification::centred);
      }
    }
  }
}

void SpectralVisualizerComponent::mouseDown(const juce::MouseEvent& e) {
  activeDragTarget = DragTarget::None;
  activeNodeIndex = -1;
  dragStartPos = e.position;

  const float w = static_cast<float>(getWidth());
  const float h = static_cast<float>(getHeight());

  // Check Transient Protection (TP) Toggle Badge click (Top-Right corner: click
  // to toggle ON/OFF in Advanced mode)
  if (isAdvancedVisible) {
    const juce::Rectangle<float> badge = getTpBadgeBounds();
    const float badgeX = badge.getX();
    const float badgeY = badge.getY();
    const float badgeW = badge.getWidth();
    const float badgeH = badge.getHeight();
    juce::Rectangle<float> badgeBounds(badgeX, badgeY, badgeW, badgeH);
    if (badgeBounds.contains(e.position)) {
      if (auto* transientParam = dynamic_cast<juce::AudioParameterBool*>(
              processor.getAPVTS().getParameter(
                  "transient_protection_enable"))) {
        *transientParam = !transientParam->get();
        transientProtectionActive = transientParam->get();
        if (!transientProtectionActive) {
          ledBrightness = 0.0f;
          transientHoldTicks = 0;
        }
        repaint();
        return;
      } else if (auto* param = processor.getAPVTS().getParameter(
                     "transient_protection_enable")) {
        float currentVal = param->getValue();
        float newVal = (currentVal > 0.5f) ? 0.0f : 1.0f;
        param->setValueNotifyingHost(newVal);
        transientProtectionActive = (newVal > 0.5f);
        if (!transientProtectionActive) {
          ledBrightness = 0.0f;
          transientHoldTicks = 0;
        }
        repaint();
        return;
      }
    }
  }

  // Check Reduction Curve nodes (if enabled and in Advanced mode)
  if (isAdvancedVisible && currentFrame.reductionCurveEnabled) {
    const auto& nodes = processor.getCurveNodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
      float nx = std::clamp(nodes[i].normX, 0.0f, 1.0f);
      float ny = SpectralVisualizerComponent::biasToY(nodes[i].biasDB, h);
      if (e.position.getDistanceFrom({nx * w, ny}) <= 12.0f) {
        activeDragTarget = DragTarget::CurveNode;
        activeNodeIndex = static_cast<int>(i);
        return;
      }
    }

    // Add a new curve node on empty click
    float clickNormX = std::clamp(e.position.x / w, 0.01f, 0.99f);
    float clickBiasDB = std::clamp(
        SpectralVisualizerComponent::yToBias(e.position.y, h), -SpectralVisualizerComponent::kCurveMaxBiasDB, SpectralVisualizerComponent::kCurveMaxBiasDB);
    processor.addCurveNode(clickNormX, clickBiasDB);
    repaint();
  }
}

void SpectralVisualizerComponent::mouseDrag(const juce::MouseEvent& e) {
  const float w = static_cast<float>(getWidth());
  const float h = static_cast<float>(getHeight());

  if (activeDragTarget == DragTarget::CurveNode && activeNodeIndex >= 0) {
    float normX = std::clamp(e.position.x / w, 0.0f, 1.0f);
    float biasDB = std::clamp(SpectralVisualizerComponent::yToBias(e.position.y, h),
                              -SpectralVisualizerComponent::kCurveMaxBiasDB, SpectralVisualizerComponent::kCurveMaxBiasDB);
    processor.updateCurveNode(activeNodeIndex, normX, biasDB);
    repaint();
  }
}

void SpectralVisualizerComponent::mouseUp(const juce::MouseEvent&) {
  activeDragTarget = DragTarget::None;
  activeNodeIndex = -1;
}

void SpectralVisualizerComponent::mouseDoubleClick(const juce::MouseEvent& e) {
  if (isAdvancedVisible && currentFrame.reductionCurveEnabled) {
    const auto& nodes = processor.getCurveNodes();
    const float w = static_cast<float>(getWidth());
    const float h = static_cast<float>(getHeight());

    for (size_t i = 1; i < nodes.size() - 1; ++i) {
      float nx = std::clamp(nodes[i].normX, 0.0f, 1.0f);
      float ny = SpectralVisualizerComponent::biasToY(nodes[i].biasDB, h);
      if (e.position.getDistanceFrom({nx * w, ny}) <= 12.0f) {
        processor.removeCurveNode(static_cast<int>(i));
        repaint();
        return;
      }
    }
  }
}

void SpectralVisualizerComponent::mouseMove(const juce::MouseEvent& e) {
  const float w = static_cast<float>(getWidth());
  juce::Rectangle<float> badgeBounds = getTpBadgeBounds();

  if (isAdvancedVisible && badgeBounds.contains(e.position)) {
    setMouseCursor(juce::MouseCursor::PointingHandCursor);
  } else {
    setMouseCursor(juce::MouseCursor::NormalCursor);
  }
}

void SpectralVisualizerComponent::mouseExit(const juce::MouseEvent&) {
  setMouseCursor(juce::MouseCursor::NormalCursor);
}

juce::String SpectralVisualizerComponent::getTooltip() {
  if (isAdvancedVisible) {
    const float w = static_cast<float>(getWidth());
    const juce::Rectangle<float> badge = getTpBadgeBounds();
    const float badgeX = badge.getX();
    const float badgeY = badge.getY();
    const float badgeW = badge.getWidth();
    const float badgeH = badge.getHeight();
    juce::Rectangle<float> badgeBounds(badgeX, badgeY, badgeW, badgeH);

    auto mousePos = getMouseXYRelative().toFloat();
    if (badgeBounds.contains(mousePos)) {
      return kTpBadgeTip;
    }
  }
  return {};
}
