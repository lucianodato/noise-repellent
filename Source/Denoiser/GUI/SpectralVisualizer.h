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

#include "../PluginProcessor.h"
#include <juce_gui_basics/juce_gui_basics.h>

class SpectralVisualizerComponent : public juce::Component,
                                    public juce::Timer,
                                    public juce::TooltipClient {
public:
  explicit SpectralVisualizerComponent(
      NoiseRepellentAudioProcessor& processorToUse);
  ~SpectralVisualizerComponent() override;

  void paint(juce::Graphics& g) override;
  void resized() override;
  void timerCallback() override;

  void mouseDown(const juce::MouseEvent& e) override;
  void mouseDrag(const juce::MouseEvent& e) override;
  void mouseUp(const juce::MouseEvent& e) override;
  void mouseDoubleClick(const juce::MouseEvent& e) override;
  void mouseMove(const juce::MouseEvent& e) override;
  void mouseExit(const juce::MouseEvent& e) override;

  juce::String getTooltip() override;

  void setAdvancedControlsVisible(bool visible) {
    if (isAdvancedVisible != visible) {
      isAdvancedVisible = visible;
      repaint();
    }
  }

  // Spectrum axis ranges (shared by mapping helpers and paint)
  static constexpr float kAxisMinDB = -100.0f;
  static constexpr float kAxisMaxDB = -20.0f;
  static constexpr float kAxisMinFreq = 20.0f;
  static constexpr float kAxisMaxFreq = 20000.0f;
  // Reduction-curve vertical mapping: +/-kCurveMaxBiasDB spans
  // kCurveHeightFrac of the height around the vertical centre
  static constexpr float kCurveMaxBiasDB = 24.0f;
  static constexpr float kCurveHeightFrac = 0.4f;
  // Transient-protection badge geometry
  static constexpr float kBadgeW = 66.0f;
  static constexpr float kBadgeH = 22.0f;
  static constexpr float kBadgeMargin = 10.0f;
  static constexpr float kBadgeTop = 10.0f;

  juce::Rectangle<float> getTpBadgeBounds() const;
  static float biasToY(float biasDB, float h);
  static float yToBias(float y, float h);

private:
  NoiseRepellentAudioProcessor& processor;
  NoiseRepellentAudioProcessor::SpectralFrame currentFrame;
  bool isAdvancedVisible = false;

  std::array<float, NoiseRepellentAudioProcessor::kFftBins> smoothedInputDB;
  std::array<float, NoiseRepellentAudioProcessor::kFftBins> smoothedOutputDB;
  bool isSmoothedInitialized = false;
  int idleTicks = 0;

  float ledBrightness = 0.0f;
  int transientHoldTicks = 0;
  bool transientProtectionActive = false;

  // Mouse Interaction
  enum class DragTarget { None, CurveNode };
  DragTarget activeDragTarget = DragTarget::None;
  int activeNodeIndex = -1;
  juce::Point<float> dragStartPos;
  float dragStartNormX = 0.0f;
  float dragStartBiasDB = 0.0f;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpectralVisualizerComponent)
};
