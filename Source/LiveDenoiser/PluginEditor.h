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

#pragma once

#include "../Shared/GUI/LookAndFeel.h"
#include "GUI/LiveSpectralVisualizer.h"
#include "PluginProcessor.h"
#include <juce_gui_basics/juce_gui_basics.h>

class NoiseRepellentLiveAudioProcessorEditor
    : public juce::AudioProcessorEditor {
public:
  explicit NoiseRepellentLiveAudioProcessorEditor(
      NoiseRepellentLiveAudioProcessor&);
  ~NoiseRepellentLiveAudioProcessorEditor() override;

  void paint(juce::Graphics&) override;
  void resized() override;

  void mouseEnter(const juce::MouseEvent&) override;
  void mouseExit(const juce::MouseEvent&) override;

private:
  NoiseRepellentLiveAudioProcessor& audioProcessor;
  NoiseRepellentLookAndFeel customLookAndFeel;

  // Header Controls
  juce::Label brandLabel;
  juce::ToggleButton btnBypass{"Bypass"};
  juce::TextButton btnDelta{"Delta"};

  // Main Reduction Control (iZotope RX Voice De-noise style)
  juce::Label lblReduction{"lblReduction", "REDUCTION"};
  juce::Slider sliderReduction{juce::Slider::LinearVertical,
                               juce::Slider::TextBoxBelow};

  // Noise floor learning (sticky toggle driving the engine tracker)
  juce::TextButton btnLearn{"Learn"};

  // Threshold bank right of display. Attack/Release/Knee are handled
  // automatically at tuned defaults (no UI, cf. RX Voice De-noise);
  // the APVTS params stay for host/session compatibility.
  juce::Slider sliderThreshold{juce::Slider::LinearVertical,
                               juce::Slider::TextBoxBelow};

  juce::Label lblThreshold{"lblThreshold", "THRESHOLD"};

  // Spectrum display
  LiveSpectralVisualizerComponent spectralVisualizer;

  // Footer Tooltip Bar
  juce::Label footerTooltipLabel;

  // Parameter Attachments
  using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
  using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;

  std::unique_ptr<ButtonAttachment> attachBypass;
  std::unique_ptr<ButtonAttachment> attachLearn;
  std::unique_ptr<SliderAttachment> attachReduction;
  std::unique_ptr<SliderAttachment> attachThreshold;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(
      NoiseRepellentLiveAudioProcessorEditor)
};
