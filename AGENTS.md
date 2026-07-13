# noise-repellent Agent Context

This file contains foundational mandates and architectural context for Gemini CLI when working on `noise-repellent`.

## Foundational Mandates

1. **RT-Safety**: This is an LV2 plugin (`hardRTCapable`). Absolutely no blocking calls in `run()`.
2. **Bypass Fidelity**: All bypass transitions must use `SignalCrossfade` to prevent pops/clicks and maintain latency alignment.
3. **Manifest Sync**: Any change to plugin parameters in `plugins/nrepellent.c` or `plugins/nrepellent-2d.c` (PortIndex, connection logic) MUST be reflected in `lv2ttl/nrepellent.ttl.in` and `lv2ttl/nrepellent-2d.ttl.in`.
   - **Parameter Cleanup**: Note that `steering_response` has been renamed to `aggressiveness`.
4. **State Persistence**: Ensure all noise profile modes (Mean, Median, Max, Min) and their respective block counts are correctly saved/restored in the LV2 state extension.

## Project Structure & Workflow

- **LV2 Wrappers**: `plugins/nrepellent.c` (1D) and `plugins/nrepellent-2d.c` (2D).
- **TTL Templates**: `.ttl.in` files located in `lv2ttl/` are used to generate the final manifest. Do not edit `.ttl` files in the build directory.
- **Subprojects**: `libspecbleach` is managed via Meson subprojects. To update the core library, modify `subprojects/libspecbleach.wrap` or configure options like `force_bundled_libspecbleach`.
- **Install & Test**: The plugin can be built and installed in release mode (linking against the local libspecbleach) using the `Install Noise Repellent` skill instructions.

## Integration Details

- **Mapping**: Note the non-linear mapping for `masking_depth` in `run()`:
  `self->parameters.masking_depth = 1.0f - powf(1.0f - (*self->masking_transparency / 100.0f), 3.0f);`
- **Latency**: Latency is reported to the host in `activate()` and through a dedicated output control port.
