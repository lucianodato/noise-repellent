---
name: meson_build
description: Best practices and commands for building, testing, and installing noise-repellent using Meson.
---

# Meson Build Skill for noise-repellent

This skill provides instructions on how to use the Meson build system for the `noise-repellent` project.

## Core Rules

1.  **Use a Dedicated Build Directory**: Always use a directory named `build` in the project root.
2.  **Avoid Redundant Build Folders**: Never create folders like `build_new`, `build2`, etc.
3.  **Reconfigure if Needed**: If the `build` directory exists but you need to change options, use:
    ```bash
    meson setup build --reconfigure [options]
    ```

## Common Workflows

### Standard Build (using system libspecbleach)
```bash
meson setup build
meson compile -C build
```

### Build with Bundled libspecbleach (Recommended for development)
```bash
meson setup build -Dforce_bundled_libspecbleach=true
meson compile -C build
```

### Build with Static libspecbleach
```bash
meson setup build -Dstatic_libspecbleach=true
meson compile -C build
```

### Custom LV2 Install Directory
```bash
# Example for macOS
meson setup build -Dlv2dir=/Library/Audio/Plug-Ins/LV2
meson compile -C build
sudo meson install -C build
```

### Debug Build with Sanitizers
```bash
meson setup build -Dbuildtype=debug -Denable_sanitizers=true
meson compile -C build
```

## Project Options

Refer to `meson_options.txt` for the full list. Key options:

| Option | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `lv2dir` | string | `''` | Install directory for LV2 bundles |
| `force_bundled_libspecbleach` | boolean | `false` | Force use of bundled libspecbleach |
| `static_libspecbleach` | boolean | `false` | Link libspecbleach statically into the plugins |
| `custom_warning_level` | integer | `2` | Compiler warning level (0-3) |
| `enable_sanitizers` | boolean | `false` | Enable sanitizers for debug builds |

## Tips
- Use `meson configure build` to see all current configuration options and their values.
- If you encounter build issues, try `ninja -C build clean` before recompiling, or as a last resort, `rm -rf build` and start over.
- RPATH is automatically handled based on whether libspecbleach is bundled, static, or system.
