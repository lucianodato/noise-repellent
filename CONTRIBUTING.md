# Contributing to Noise Repellent

We welcome your contributions! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features

## Development Workflow

We use **CMake** and **JUCE** for building.

1. **Clone the repo:**
   ```bash
   git clone https://github.com/lucianodato/noise-repellent.git
   cd noise-repellent
   ```

2. **Setup build:**
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Debug
   ```

3. **Compile:**
   ```bash
   cmake --build build --config Debug -j4
   ```

4. **Format Code:**
   We use `clang-format` according to the project's `.clang-format` rules. Please format modified files before submitting:
   ```bash
   clang-format -i Source/PluginProcessor.cpp Source/PluginProcessor.h
   ```

## Pull Requests

1. Fork the repo and create your branch from `master` or `main`.
2. Ensure existing functionality builds cleanly across target plugin formats.
3. If you've added features or modified UI components, verify behavior in test hosts.
4. Make sure your code complies with our coding guidelines.
5. Open a pull request!

## License

By contributing, you agree that your contributions will be licensed under the GNU General Public License v3.0 (GPL-3.0-or-later).
