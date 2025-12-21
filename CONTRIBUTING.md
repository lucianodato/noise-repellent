# Contributing to Noise Repellent

We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features

## Development Workflow

We use **Meson** and **Ninja** for building.

1.  **Clone the repo and subprojects:**
    ```bash
    git clone --recursive https://github.com/lucianodato/noise-repellent.git
    cd noise-repellent
    ```

2.  **Setup build:**
    ```bash
    meson setup build --buildtype=debug -Denable_sanitizers=true
    ```

3.  **Compile and Test:**
    ```bash
    meson compile -C build
    ```

4.  **Format Code:**
    We use `clang-format`. Please format your code before submitting:
    ```bash
    ninja -C build format
    ```

## Pull Requests

1.  Fork the repo and create your branch from `master`.
2.  If you've added code that should be tested, add tests.
3.  If you've changed APIs, update the documentation.
4.  Ensure the test suite passes.
5.  Make sure your code lints.
6.  Issue that pull request!

## License

By contributing, you agree that your contributions will be licensed under its LGPL-3.0 License.
