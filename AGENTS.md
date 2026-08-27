# noise-repellent Agent Context

This file contains foundational mandates and architectural context for AI agents working on `noise-repellent`.

## Foundational Mandates

1. **RT-Safety**: `processBlock()` executes in the real-time audio thread. Strictly NO dynamic allocations (`malloc`, `new`), NO blocking primitives (locks, mutexes, condition variables), and NO file/console I/O.
2. **Bypass Fidelity**: All bypass transitions must use native `juce::dsp::DryWetMixer<float>` to prevent pops/clicks and maintain latency alignment (`setWetLatency`) with host delay compensation.
3. **APVTS & Parameter Management**: Parameter changes must go through `juce::AudioProcessorValueTreeState`. Parameter IDs defined in `PluginProcessor` layout must match key strings used in state serialization and GUI attachments.
4. **License**: Code in `Source/` is licensed under GPL-3.0-or-later. Include standard license headers on all new source files.
5. **UX & Invariant Preservation**: All changes must preserve the contracts defined in [docs/UX_INVARIANTS.md](docs/UX_INVARIANTS.md). Always verify with `test_ux_invariants` (`-DENABLE_PLUGIN_TESTS=ON`).

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

## Developer & AI Agent Guidelines for NoiseRepellent

This document defines mandatory rules for maintaining, building, and releasing NoiseRepellent. Any AI agent or developer making changes to `CMakeLists.txt`, build workflows, or project metadata MUST follow these constraints.

---

## 1. Versioning & Release Workflow

* **Single Source of Truth**: `CMakeLists.txt` defines the project version via `project(NoiseRepellent VERSION X.Y.Z LANGUAGES C CXX)`.
* **Release Process**:
  1. Bump `VERSION` in `CMakeLists.txt` (e.g., `0.3.1`).
  2. Commit and push the version bump to `master`.
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

## 3. CMake & Dependency Linking Rules

To balance portable DAW binary compatibility for release artifacts with downstream Linux distribution packaging standards:

1. **Dual-Mode Dependency Linking**:
   * **Default Portable Release Mode (`USE_SYSTEM_*=OFF`)**:
     * Third-party libraries (`FreeType`, `libspecbleach`) **must be statically embedded** by default into plugin targets to ensure zero host DAW crashes or symbol conflicts (e.g., in Ardour, Bitwig, or REAPER).
     * On Linux, static GCC runtimes (`-static-libgcc -static-libstdc++`) must be linked for release builds.
     * On Windows, static MSVC runtimes (`/MT`) must be configured (`set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")`).
   * **Downstream Packaging Mode (`USE_SYSTEM_*=ON`)**:
     * Provide CMake build toggles (`USE_SYSTEM_FREETYPE`, `USE_SYSTEM_SPECBLEACH`, `USE_SYSTEM_JUCE`) defaulting to `OFF`.
     * When set to `ON` by Linux distro packagers, CMake uses `find_package()` / `pkg_check_modules()` to dynamically link against system shared libraries (`.so`).
2. **Package Metadata Requirements**:
   * Standard prebuilt release package configurations (`.deb` / CPack) declare core desktop system dependencies (`libasound2`, `libgl1`, `libx11-6`, `libxext6`, `libxcursor1`, `libxinerama1`, `libxrandr2`). Distro-packaged builds managed by package maintainers will handle full dynamic dependency trees via `dpkg-shlibdeps`.
3. **Position Independent Code**:
   * `set(CMAKE_POSITION_INDEPENDENT_CODE ON)` must remain enabled globally so static `.a` libraries compile with `-fPIC` for shared object (`.so` / `.vst3` / `.dylib`) plugin outputs.
4. **Local Submodule Development**:
   * The CMake build must automatically detect local copies of `libspecbleach` at `../libspecbleach` or `subprojects/libspecbleach` before falling back to remote fetching.
   * Internal subproject source properties (e.g., `nlm_filter_avx.c`) must attach directly to the `libspecbleach` target so AVX symbols resolve cleanly across x86_64 architectures.

---

## 4. CI/CD & Binary Verification Requirements

The GitHub Actions workflow enforces strict post-build binary verification steps across all OS targets:

* **Linux**: For standard release builds, runs `ldd` against a strict whitelist of core OS/desktop runtimes (`libc`, `libm`, `libpthread`, `libdl`, `librt`, `ld-linux`, `libasound`, `libGL`, `libX11`, `libXext`, `libXcursor`, `libXinerama`, `libXrandr`, `libXrender`, `libXi`, `libgomp`). Fails if unexpected dynamic libraries (`libfreetype.so`, etc.) or dynamic C++ runtimes (`libstdc++.so`, `libgcc_s.so`) leak into the prebuilt release binaries.
* **macOS**: Uses `lipo` to verify fat universal binaries (`arm64` + `x86_64`) and `otool -L` to ensure zero Homebrew or non-system dynamic library leaks.
* **Windows**: Parses PE import tables via `objdump -p` to guarantee only standard system DLLs (`KERNEL32.dll`, `VCRUNTIME140.dll`, `ucrtbase.dll`, `d2d1.dll`, etc.) are imported.

