# Changelog

All notable changes to this project will be documented in this file.

## [0.3.2] - Unreleased

### Added
- **Improve Adaptive/Manual profile learn UX and UI**: Adaptive noise learn can work as a standalone mode or on top of a manual profile.

### Improved & Refactored

### Fixed

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
- Fixed bypass issue for hosts using in-place buffers (e.g., Ardour) (Thanks @jmaibaum).
- Fixed soft bypass latency compensation (Issue #124) and stereo crossfade independence.

### Added
- New control parameters from development branch:
    - **Noise Scaling Type**: Choose between different reduction algorithms.
    - **Post-filter Threshold**: Fine-tune the reduction threshold.
    - **Residual Whitening**: Added to 2D denoiser (replacing Over-subtraction), matching standard plugin behavior.
