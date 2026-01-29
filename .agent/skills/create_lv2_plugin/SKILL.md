---
name: create_lv2_plugin
description: Instructions for creating a new LV2 plugin in noise-repellent using libspecbleach.
---

# Create LV2 Plugin Skill

This skill guides you through wrapping a `libspecbleach` processor as an LV2 plugin in `noise-repellent`.

## 1. Plugin Implementation
Location: `plugins/nrepellent-<name>.c`

1.  **Struct**: Define `NRepellent<Name>` struct.
    *   Include `SpectralBleachHandle` from `libspecbleach`.
    *   Include `float*` pointers for ports.
2.  **Ports**: Define `PortIndex` enum (Audio In/Out, Controls).
3.  **Instantiate**:
    *   Map `libspecbleach` initialization.
    *   Handle sample rate and frame size.
4.  **Connect Port**: Map LV2 buffers to struct pointers.
5.  **Run**:
    *   Read control port values.
    *   Map to `SpectralBleachParameters`.
    *   Call `specbleach_<name>_process`.
    *   **Latency**: If the processor has latency, report it via LV2 meta (if dynamic) or ensure it's documented.
6.  **Cleanup**: Free the library handle.

## 2. TTL Definitions
Location: `lv2ttl/`

1.  Create `nrepellent-<name>.ttl.in`.
    *   Define the plugin URI (e.g., `http://github.com/lucianodato/noise-repellent#<name>`).
    *   Define all ports (inputs, outputs, controls) with correct indices.
2.  Update `manifest.ttl.in` to include the new plugin URI and binary.

## 3. Build System
Location: `meson.build`

1.  Add `nrepellent-<name>` target using `shared_module`.
2.  Add `nrepellent-<name>.ttl` configuration step.

## 4. Verification
1.  **Build**: `meson setup build --buildtype=release` then `meson compile -C build`.
    *   **CRITICAL**: Always use Release mode (`-O3`) for verification. Debug builds will cause xruns with heavy DSP.
2.  **Verify**: Ensure no link errors with `libspecbleach`.
