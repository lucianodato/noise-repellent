# Noise Repellent

[![CI Build](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/lucianodato/noise-repellent?utm_source=oss&utm_medium=github&utm_campaign=lucianodato%2Fnoise-repellent&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

A multi-format audio plugin (VST3, AU, LV2) for real-time spectral noise reduction, built with [JUCE](https://juce.com/) and powered by the [libspecbleach](https://github.com/lucianodato/libspecbleach) DSP engine.

## Screenshots
![Advanced Controls](<Images/Screenshot 1.png>)
![Basic Controls](<Images/Screenshot 2.png>)

## Features

### Advanced Denoising Algorithms
* **2D Non-Local Means (NLM)**: Uses spectral-temporal pattern matching to suppress musical noise while preserving high-frequency detail and textures.
* **Manual Profiling**: Classic noise reduction using a user-captured noise profile from a silent section.
* **Adaptive Estimation**: Real-time noise floor tracking with multiple algorithms.
    * **Hybrid Operation**: Works on top of manual profiles to refine captured snapshots in real-time.
    * **SPP-MMSE**: Robust, unbiased estimation for complex noise environments (best for voice).
    * **Brandt (Trimmed Mean)**: Efficient estimation for steady-state broadband noise.
    * **Martin Minimum Statistics**: Reliable tracking for slowly varying noise.

### Precision Controls
* **Tonal separation**: Independent reduction of harmonic content and tonal noise (hum, resonance).
* **Intelligent Steering**: Adjustable aggressiveness to balance between different noise profile statistics (Mean, Median, Max).
* **Adaptive Whitening**: Reshapes the residual noise floor to prevent coloring and artifacts.
* **Masking Transparency (Veto)**: Protects transients and delicate details by balancing reduction against signal energy.
* **Adjustable Smoothing**: Fine-tune the balance between artifact suppression and transient clarity.
    * **Standard Denoising**: Uses temporal smoothing for stable noise reduction.
    * **2D Denoising**: Uses NLM smoothing for pattern-based artifact removal.

### Workflow & Integration
* **Interactive Spectral Visualizer**: Real-time FFT visualization of input, noise floor profile, and processed output spectrums with detected tonal peak markers.
* **Residual Listening**: Hear exactly what is being filtered out to fine-tune your settings.
* **Soft Bypass**: Seamless, click-free A/B testing with cross-faded bypass and latency compensation.
* **Full State Saving**: Noise profiles and all APVTS parameters are saved with the host session.

### Multi-Format Compatibility
* **Formats**: VST3, AU, and LV2.
* **Platforms**: Optimized for Linux, macOS, and Windows.
* **Stereo Support**: Full multi-channel processing ready for modern stereo production workflows.

## Installation

Pre-built installers and packages for Linux, macOS, and Windows are available on the [GitHub Releases](https://github.com/lucianodato/noise-repellent/releases) page.

### 🪟 Windows
1. Download `NoiseRepellent-Win64-Installer.exe` from the latest release.
2. Run the installer to automatically place the VST3 and LV2 formats into your standard System plugin directories (`C:\Program Files\Common Files\`).

---

### 🍎 macOS (Universal: Apple Silicon & Intel)
1. Download `NoiseRepellent-macOS-Universal.pkg`.
2. **Right-Click** (or **Ctrl + Click**) the `.pkg` file and choose **Open** from the menu to bypass the initial Gatekeeper prompt.
3. Follow the installation wizard.

> [!NOTE]
> The macOS installer automatically strips Gatekeeper quarantine flags during installation, so your installed VST3, AU, and LV2 formats will load immediately inside your DAW without requiring any Terminal commands!

---

### 🐧 Linux (x86_64)

#### Debian / Ubuntu / Linux Mint / Pop!_OS
1. Download `noise-repellent_0.3.1_amd64.deb`.
2. Install via software center or terminal:
   ```bash
   sudo apt install ./noise-repellent_0.3.1_amd64.deb

#### Other Linux Distributions
1. Download `noise-repellent-linux-x86_64.tar.gz`
2. tar -xvf noise-repellent-linux-x86_64.tar.gz
3. cd noise-repellent-linux-x86_64
4. ./install.sh

### From Source

**Requirements:**
- CMake >= 3.22
- C11 / C++17 Compiler (GCC 10+, Clang 12+, or MSVC 2019+)
- Ninja or Make
- Linux GUI build dependencies (Linux only):
  - `libasound2-dev`, `libgl1-mesa-dev`, `libx11-dev`, `libxcomposite-dev`, `libxcursor-dev`, `libxext-dev`, `libxinerama-dev`, `libxrandr-dev`, `libfontconfig1-dev`, `libfreetype6-dev`

---

#### 1. Standard Portable Build (Self-Contained)
By default, CMake fetches and statically embeds dependencies (`FFTW3`, `FreeType`, `libspecbleach`) to produce standalone, portable release binaries that prevent symbol conflicts across different DAWs and Linux distributions.

```bash
# Clone the repository including submodules
git clone --recurse-submodules https://github.com/lucianodato/noise-repellent.git
cd noise-repellent

# Configure build with static bundled dependencies
cmake -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release

# Compile all plugin formats (VST3, AU, LV2)
cmake --build build --config Release --parallel

# Install plugins to standard system plugin locations
sudo cmake --install build
```

#### 2. Downstream Linux Packaging Build (Shared System Libraries)
Linux distribution packagers (Arch, Debian, Fedora, etc.) can configure the build to dynamically link against host system libraries (`libfftw3f`, `libfreetype`, `libspecbleach`) using CMake build options.

> **Note:** Requires system development packages installed on the host (e.g., `libfftw3-dev`, `libfreetype6-dev`, `libspecbleach-dev`).

```bash
# Configure build using system-installed shared libraries
cmake -B build -G "Ninja" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr \
  -DUSE_SYSTEM_FFTW=ON \
  -DUSE_SYSTEM_FREETYPE=ON \
  -DUSE_SYSTEM_SPECBLEACH=ON

# Compile all plugin formats
cmake --build build --config Release --parallel

# Install to packaging staging area or system
cmake --install build
```

The compiled plugin targets (VST3, AU, LV2) will be generated in `build/NoiseRepellent_artefacts/Release/` and installed to system plugin locations (`/usr/lib/vst3`, `/usr/lib/lv2`, `/Library/Audio/Plug-Ins/`, etc.).

#### Custom Installation Paths
You can override the installation destination for specific formats using `-Doption=VALUE`:

| Option | Default | Description |
| :--- | :--- | :--- |
| `INSTALL_VST3_DIR` | OS Standard (`${CMAKE_INSTALL_LIBDIR}/vst3` on Linux, `/Library/Audio/Plug-Ins/VST3` on macOS) | Destination path for VST3 plugin bundle |
| `INSTALL_AU_DIR` | OS Standard (`/Library/Audio/Plug-Ins/Components` on macOS) | Destination path for AU component bundle |
| `INSTALL_LV2_DIR` | OS Standard (`${CMAKE_INSTALL_LIBDIR}/lv2` on Linux, `/Library/Audio/Plug-Ins/LV2` on macOS) | Destination path for LV2 plugin bundle |

> [!IMPORTANT]
> Noise Repellent relies heavily on real-time FFT processing and AVX vectorization routines. Always configure CMake with `-DCMAKE_BUILD_TYPE=Release` to prevent CPU spikes and audio buffer xruns in your DAW.


## Usage

Please refer to the [Project Wiki](https://github.com/lucianodato/noise-repellent/wiki) for detailed usage instructions.

## Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our development workflow and standards.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0-or-later) - see the [LICENSE](LICENSE) file for details.
