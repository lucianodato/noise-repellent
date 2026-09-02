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

  // Adaptive two-part button (toggle + method dropdown, like the Denoiser)
  juce::TextButton btnAdaptiveNoise{"Adaptive"};
  juce::TextButton btnAdaptiveArrow{juce::CharPointer_UTF8("\xe2\x96\xbc")};

  // Noise floor learning (GUI-only state, no automation)
  juce::TextButton btnLearn{"Learn"};

  // Advanced Controls Panel (Collapsible)
  juce::TextButton btnAdvancedToggle{"ADVANCED"};
  juce::GroupComponent groupAdvanced{"groupAdvanced", "ADVANCED CONTROLS"};
  // Invisible state holder for the adaptive_method attachment
  juce::ComboBox comboMethod;
  juce::Slider sliderSmoothing{juce::Slider::RotaryHorizontalVerticalDrag,
                               juce::Slider::NoTextBox};
  juce::Slider sliderSuppression{juce::Slider::RotaryHorizontalVerticalDrag,
                                 juce::Slider::NoTextBox};

  juce::Label lblSmoothing{"lblSmoothing", "SMOOTHING"};
  juce::Label lblSuppression{"lblSuppression", "AGGRESSIVENESS"};

  // Spectrum display
  LiveSpectralVisualizerComponent spectralVisualizer;

  // Footer Tooltip Bar
  juce::Label footerTooltipLabel;

  // Parameter Attachments
  using ButtonAttachment = juce::AudioProcessorValueTreeState::ButtonAttachment;
  using SliderAttachment = juce::AudioProcessorValueTreeState::SliderAttachment;
  using ComboBoxAttachment =
      juce::AudioProcessorValueTreeState::ComboBoxAttachment;

  std::unique_ptr<ButtonAttachment> attachBypass;
  std::unique_ptr<ButtonAttachment> attachAdaptive;
  std::unique_ptr<SliderAttachment> attachReduction;
  std::unique_ptr<ComboBoxAttachment> attachMethod;
  std::unique_ptr<SliderAttachment> attachSmoothing;
  std::unique_ptr<SliderAttachment> attachSuppression;

  bool isAdvancedVisible = false;

  void updateLayout();

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(
      NoiseRepellentLiveAudioProcessorEditor)
};
