# noise-repellent Agent Context

This file contains foundational mandates and architectural context for AI agents working on `noise-repellent`.

## Foundational Mandates

1. **RT-Safety**: `processBlock()` executes in the real-time audio thread. Strictly NO dynamic allocations (`malloc`, `new`), NO blocking primitives (locks, mutexes, condition variables), and NO file/console I/O.
2. **Bypass Fidelity**: All bypass transitions must use native `juce::dsp::DryWetMixer<float>` to prevent pops/clicks and maintain latency alignment (`setWetLatency`) with host delay compensation.
3. **APVTS & Parameter Management**: Parameter changes must go through `juce::AudioProcessorValueTreeState`. Parameter IDs defined in `PluginProcessor` layout must match key strings used in state serialization and GUI attachments.
4. **License**: Code in `Source/` is licensed under GPL-3.0-or-later. Include standard license headers on all new source files.

## Project Structure & Workflow

- **Plugin Core**: `Source/PluginProcessor.h` & `Source/PluginProcessor.cpp` manage audio processing, APVTS parameters, libspecbleach handles (`specbleach1D`, `specbleach2D`), `juce::dsp::DryWetMixer`, and FFT visualization fifo.
- **GUI Engine**: `Source/PluginEditor.h` & `Source/PluginEditor.cpp` assemble the layout.
- **GUI Components**: `Source/GUI/LookAndFeel.h` / `.cpp` for visual styling and `Source/GUI/SpectralVisualizer.h` / `.cpp` for FFT visualization.
- **Build & Dependencies**: Uses CMake 3.22+. JUCE 8 is fetched via `FetchContent`. `libspecbleach` is linked as a static library via `add_subdirectory`.

## Build Commands

```bash
# Debug build
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug -j4

# Release build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j4
```

## Integration Details

- **Bypass**: Soft bypass is handled via native `juce::dsp::DryWetMixer<float>`, configured with `setWetLatency(latency)` for dry delay compensation and `setWetMixProportion()` for smooth crossfading. `processBlockBypassed()` is also implemented.
- **FFT Visualization**: Input and output channels feed into `juce::dsp::FFT` order 12 (4096 points) and push magnitude spectra to an SPSC lock-free `juce::AbstractFifo` for smooth rendering in `SpectralVisualizer`.
- **Latency**: Plugin latency is updated dynamically based on algorithm mode (1D STFT vs 2D NLM) via `setLatencySamples()`.
