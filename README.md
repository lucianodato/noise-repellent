# Noise Repellent

A suite of LV2 plugins for real-time spectral noise reduction, built on the [libspecbleach](https://github.com/lucianodato/libspecbleach) library.

[![CI Build](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml)

## Features

* **Manual capture noise reduction**: Typical noise reduction where you capture a noise profile.
* **Adaptive noise reduction**: Automatic noise suppression, optimized for voice (low latency).
* **Adjustable parameters**: Reduction amount, smoothing, whitening, and transient protection.
* **Residual listening**: Hear exactly what is being removed.
* **Soft bypass**: Cross-faded bypass to avoid clicks.
* **State saving**: Noise profiles are saved with your host session.

## Screenshots

![Noise Repellent in Ardour](images/Noise-repellent-Ardour.png)
![Noise Repellent in Reaper](images/Noise-repellent-Reaper.png)
![Adaptive Noise Repellent in Reaper](images/Adaptive-Reaper.png)

## Installation

### From Binaries
Binaries for Linux, macOS, and Windows are provided in the [GitHub Releases](https://github.com/lucianodato/noise-repellent/releases) page. Extract the folder to your LV2 plugins directory:
- **Linux**: `~/.lv2/` or `/usr/lib/lv2/`
- **macOS**: `~/Library/Audio/Plug-Ins/LV2/` or `/Library/Audio/Plug-Ins/LV2/`
- **Windows**: `%COMMONPROGRAMFILES%\LV2\`

### From Source

**Requirements:**
- Meson build system >= 0.60.0 & Ninja
- C Compiler (GCC/Clang)
- LV2 development headers
- `pkg-config`

**Build:**

```bash
git clone https://github.com/lucianodato/noise-repellent.git
cd noise-repellent

# Configure build
meson setup build --buildtype=release

# Compile
meson compile -C build

# Install (sudo may be required)
meson install -C build
```

**Build Options:**

You can configure the build options using `-Doption=value`:

- `custom_warning_level`: 0-3 (default: 2). Controls compiler warning verbosity.
- `treat_warnings_as_errors`: Treat compiler warnings as errors (default: false).
- `enable_sanitizers`: Enable sanitizers for debug builds (default: false).
- `sanitize_address`: Enable address sanitizer (only if enable_sanitizers is true) (default: true).
- `sanitize_undefined`: Enable undefined behavior sanitizer (only if enable_sanitizers is true) (default: true).
- `lv2dir`: Install directory for LV2 bundles (absolute path or relative to prefix) (default: '').
- `force_bundled_libspecbleach`: Force use of bundled libspecbleach instead of system version (default: false). Enable this to ensure API compatibility when building from source.

Example:
```bash
meson setup build -Dbuildtype=debug
```

## Usage

Please refer to the [Project Wiki](https://github.com/lucianodato/noise-repellent/wiki) for detailed usage instructions.

## Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## License

This project is licensed under the LGPL-3.0 License - see the [LICENSE](LICENSE) file for details.
