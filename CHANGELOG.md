# Changelog

All notable changes to this project will be documented in this file.

## [0.4.0] - 2026-09-07

### Added
- **Low-latency mode**: New non-automatable `low_latency` option (Options menu) for live scenarios: fixed 512-sample frame (~10.7 ms at 48 kHz / ~11.6 ms at 44.1 kHz), causal 1D-only engine, smoothing selector locked to Standard, frame-size menu disabled, and reduction/threshold links forced on. Entering it installs control defaults (smoothing floor 30, aggressiveness 1, masking 0) and starts clean (profile dropped, Learn auto-stopped). The smoothing slider's top half maps to ~130 ms max release (long releases smear coarse LF bins into pad artifacts). Toggling suspends, rebuilds from a clean slate, and re-reports PDC.
- **Third algorithm mode**: Patch-Based + Refinement (post-NLM DFTT) alongside Standard (1D) and Patch-Based (2D NLM); library-owned gapless transitions with allocation-free internal crossfade.
- **Stepped STFT frame sizes**: Options menu offers 23 / 32 / 46 / 64 / 93 ms frames. Switching suspends, rebuilds from a clean slate (profile dropped, Learn auto-stopped), and re-reports PDC; session state restores are exempt.
- **Unlinkable tonal threshold offset**: `link_threshold_offset` toggle allows independent tonal vs broadband threshold control, with synthetic tonal-mask synthesis keeping the UI responsive before the DSP pipeline runs.
- **Threshold offset & custom reduction curve**: User-controllable offsets and curve mapping.
- **Transient protection toggle & quality selection**: Transient-protection switch on the visualizer LED with latency-compensated transitions.
- **Live profile rendering**: Noise profile renders in real time while learning; engine queried for tonal peaks to update aggressiveness thresholds.
- **Adaptive/manual learn UX**: Adaptive noise learn works standalone or on top of a manual profile; refined profile UX with compact advanced panel and simplified default experience.
- **Offline-render detection**: UI overlay indicates offline rendering state.

### Improved & Refactored
- **Adopted redesigned libspecbleach C API**: Type-safe handles, extras orchestration layer (stereo groups, transitions), engine sync moved off the audio thread, unified denoiser with library-owned mode switching.
- **PFFFT / OpenMP removal**: FFTW3 dependency gone (vendored PFFFT in the library); Windows builds no longer ship non-redistributable OpenMP runtimes, and libspecbleach embeds cleanly in DLLs.
- **Bypass fidelity**: DSP skipped when bypassed via native DryWetMixer with latency-compensated crossfades; engine keeps running so toggles never time-travel.
- **Gapless engine switch**: Deferred PDC reporting with blocking overlay during rebuilds.
- **Controls**: Suppression parameter removed, aggressiveness slider moved; layout updates on algorithm mode change with improved link/smoothing coordination.

### Fixed
- Frame-size description clarity and bypass alignment issues.

**Note**: Release binaries are built by CI, which also produces the macOS universal (arm64 + x86_64) artifacts via lipo; local `cmake -B build` yields a single-arch dev build.

## [0.3.1] - 2026-08-08

### Changed
- **Major Architecture Migration**: Migrated plugin framework from legacy LV2 C implementation to modern C++ framework using JUCE 8.
- **Multi-Format Support**: Now builds natively as VST3, AU, and LV2 audio plugins for macOS, Linux, and Windows.
- **Build System**: Migrated build system from Meson/Ninja to CMake 3.22+.
- **License Update**: License updated to GNU General Public License v3.0 (GPL-3.0-or-later) to comply with JUCE open-source licensing.
- **Modern DSP Wrapper**: Refactored internal DSP wrappers (e.g., `SignalCrossfade`) from C struct helpers to modern C++ classes under the `noise_repellent` namespace.

### Added
- **Interactive Spectral Visualizer**: Custom JUCE GUI component featuring real-time FFT spectrum display (input, noise floor profile, and output) with overlay markers for detected tonal peaks.
- **APVTS State Persistence**: Complete parameter and noise profile binary state persistence using JUCE AudioProcessorValueTreeState.

## [0.3.0] - 2026-01-28

### Added
- Modern build system using Meson and Ninja.
- Configurable default frame size via `default_frame_size_ms` build option.
- Code formatting target using `clang-format`.
- CI/CD workflow for automated building and testing on Linux, macOS, and Windows.
- Detailed README and Contributing guide.

### Changed
- Updated internal API to match modern `libspecbleach` (added frame size parameter).
- Default frame size is now 40ms to support low-latency requirements.
- Updated compiler flags to strict C17 standard usage.

### Fixed
- **Critical**: Fixed audio alignment/latency issue in 2D denoiser (caused by NLM look-ahead mismatch).
- Fixed strict prototype warnings in C headers.
- Fixed stereo state restoration logic (Thanks @orivej).
- Fixed bypass issue for hosts using in-place buffers (Thanks @jmaibaum).
- Fixed soft bypass latency compensation (Issue #124) and stereo crossfade independence.

### Added
- New control parameters from development branch:
    - **Noise Scaling Type**: Choose between different reduction algorithms.
    - **Post-filter Threshold**: Fine-tune the reduction threshold.
    - **Residual Whitening**: Added to 2D denoiser (replacing Over-subtraction), matching standard plugin behavior.
