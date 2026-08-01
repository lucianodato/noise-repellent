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

#include <juce_gui_basics/juce_gui_basics.h>
#include "../PluginProcessor.h"

class SpectralVisualizerComponent : public juce::Component,
                                    public juce::Timer
{
public:
    explicit SpectralVisualizerComponent(NoiseRepellentAudioProcessor& processorToUse);
    ~SpectralVisualizerComponent() override;

    void paint(juce::Graphics& g) override;
    void resized() override;
    void timerCallback() override;

private:
    NoiseRepellentAudioProcessor& processor;
    NoiseRepellentAudioProcessor::SpectralFrame currentFrame;

    std::array<float, NoiseRepellentAudioProcessor::kFftBins> smoothedInputDB;
    std::array<float, NoiseRepellentAudioProcessor::kFftBins> smoothedOutputDB;
    bool isSmoothedInitialized = false;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpectralVisualizerComponent)
};
