# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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
- Fixed strict prototype warnings in C headers.
- Fixed stereo state restoration logic (Thanks @orivej).
- Fixed bypass issue for hosts using in-place buffers (e.g., Ardour) (Thanks @jmaibaum).
- Fixed soft bypass latency compensation (Issue #124) and stereo crossfade independence.

### Added
- New control parameters from development branch:
    - **Noise Scaling Type**: Choose between different reduction algorithms.
    - **Post-filter Threshold**: Fine-tune the reduction threshold.
    - **Residual Whitening**: Added to the standard plugin (previously only in adaptive).
