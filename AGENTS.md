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

# Developer & AI Agent Guidelines for NoiseRepellent

This document defines mandatory rules for maintaining, building, and releasing NoiseRepellent. Any AI agent or developer making changes to `CMakeLists.txt`, build workflows, or project metadata MUST follow these constraints.

---

## 1. Versioning & Release Workflow

* **Single Source of Truth**: `CMakeLists.txt` defines the project version via `project(NoiseRepellent VERSION X.Y.Z LANGUAGES C CXX)`.
* **Release Process**:
  1. Bump `VERSION` in `CMakeLists.txt` (e.g., `0.3.1`).
  2. Commit and push the version bump to `main`.
  3. Create and push a Git tag formatted as `vX.Y.Z` (e.g., `v0.3.1`).
* **CI Version Resolution**: The GitHub Actions workflow dynamically extracts the version from `git describe` and populates `PLUGIN_VERSION`. Do not hardcode version fallback strings inside `.github/workflows/`.

---

## 2. LV2 URI Guidelines

* **Current URI**: `https://github.com/lucianodato/noise-repellent#v0`
* **Fragment Anchor (`#v0`)**: Uses the URL fragment trick so the URI is both a valid, clickable GitHub link in browsers AND a unique string ID for DAWs.
* **Compatibility Rule**: 
  * **DO NOT** change or bump the LV2 URI during minor or patch updates (e.g., `0` $\rightarrow$ `1`). Keeping the URI static preserves DAW session state compatibility.
  * **ONLY** change the URI if a major breaking redesign introduces incompatible parameter structures or state formats.

---

## 3. CMake & Static Linking Rules

To ensure prebuilt binaries work seamlessly across all host DAWs and Linux distributions without symbol missing errors (like `FT_Get_Paint`):

1. **Zero External Dynamic Dependencies**:
   * All third-party libraries (`FFTW3`, `FreeType`, `libspecbleach`) **must be statically embedded** into the plugin targets via CMake `FetchContent`.
   * On Linux, static GCC runtimes (`-static-libgcc -static-libstdc++`) must be linked.
   * On Windows, static MSVC runtimes (`/MT`) must be used (`set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")`).
2. **Order of Evaluation**:
   * Declare static `freetype` **BEFORE** `FetchContent_MakeAvailable(JUCE)` on Linux to prevent JUCE from auto-linking system `libfreetype.so`.
3. **Position Independent Code**:
   * `set(CMAKE_POSITION_INDEPENDENT_CODE ON)` must remain enabled globally so static `.a` libraries compile with `-fPIC` for dynamic `.so` plugin output.
4. **Local Submodule Development**:
   * The CMake build must automatically detect local copies of `libspecbleach` at `../libspecbleach` or `subprojects/libspecbleach` before falling back to remote fetching.
   * Internal subproject source properties (e.g., `nlm_filter_avx.c`) must attach directly to `libspecbleach` target so AVX symbols resolve cleanly across x86_64 architectures.

---

## 4. CI/CD & Binary Verification Requirements

The GitHub Actions workflow enforces strict post-build binary verification steps across all OS targets:

* **Linux**: Runs `ldd` and `objdump` to enforce a strict whitelist (`libc`, `libm`, `libpthread`, `libdl`, `librt`, `ld-linux`). Fails if unexpected libraries (`libfreetype.so`, `libfftw3.so`) or undefined symbols (`FT_`, `fftw_`) leak into the binary.
* **macOS**: Uses `lipo` to verify fat universal binaries (`arm64` + `x86_64`) and `otool -L` to ensure zero Homebrew/third-party dylib leaks.
* **Windows**: Parses PE import tables via `objdump -p` to guarantee only standard system DLLs (`KERNEL32.dll`, `VCRUNTIME140.dll`, etc.) are imported.
