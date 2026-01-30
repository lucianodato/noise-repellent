---
name: development_guidelines
description: "LV2 plugin development, maintenance, and project architecture guidelines for noise-repellent"
---

# noise-repellent Development Guidelines

## LV2 Plugin Development Guidelines

**CRITICAL**: Always rebuild and test after modifying plugin code or TTL files.

### Plugin Architecture

#### Plugin Types
- **nrepellent**: Standard noise reduction plugin (manual noise profile)
- **nrepellent-adaptive**: Adaptive noise reduction plugin (automatic noise learning)

### Common Port Operations
See the [Port Management Workflow](file://../../workflows/port_management.md) for details on adding/removing ports.

### TTL File Management

#### File Structure
- **`.ttl.in` files**: Templates in `lv2ttl/` directory
- **Generated `.ttl` files**: Built in `build/` directory
- **Never edit generated `.ttl` files directly**

### Plugin Registration
Each plugin needs:
- URI definition (e.g., `NOISEREPELLENT_ADAPTIVE_URI`)
- Descriptor struct (`LV2_Descriptor`)
- `lv2_descriptor()` function entry

## Build System

### Meson Configuration
```bash
# Configure build
meson setup build --wrap-mode=forcefallback

# Build project
meson compile -C build

# Clean rebuild
rm -rf build && meson setup build --wrap-mode=forcefallback
```

## LV2 Plugin Development Patterns

### Plugin Structure Requirements

#### Required Functions
- `instantiate()`: Create plugin instance
- `connect_port()`: Connect control/audio ports
- `activate()`: Prepare for processing
- `run()`: Process audio (required)
- `deactivate()`: Optional cleanup
- `cleanup()`: Destroy instance

### Memory Management
- Use `calloc()` for zero-initialized memory
- Check allocation success
- Store pointers in plugin struct
- Free all allocated memory in `cleanup()`

### Audio Processing
- Support both separate and in-place buffers
- Copy input buffers when in-place processing needed
- Use crossfade for smooth parameter transitions
- Store sample rate from `instantiate()`

### Error Handling
- Check all required features in `instantiate()`
- Validate allocation success
- Log errors with `lv2_log_error()`

### Plugin Variants
- Separate descriptor functions for each variant
- Share common processing code
- Define unique URIs for each variant

---

## LV2 Plugin Safety and Reliability Rules

### Control Port Safety
**CRITICAL**: Always protect against unconnected control ports in `run()` functions.

**Pattern**:
```c
self->parameters.param = self->param_port ? *self->param_port : DEFAULT_VALUE;
```

### Audio Port Validation
**CRITICAL**: Validate audio ports before processing.

### Latency Compensation Handling
**CRITICAL**: Handle initial latency periods in crossfade and delay compensation logic.

### Processing State Management
**IMPORTANT**: Always call core processing functions regardless of UI state.

### Parameter Validation
**RECOMMENDED**: Validate parameter ranges before passing to processing functions.

---

## Quality Assurance

### Build Verification
- [ ] All plugins compile without warnings
- [ ] TTL files generate correctly
- [ ] Port indices are synchronized
- [ ] Plugin binaries load in test host

### Code Review Checklist
- [ ] Port management follows guidelines
- [ ] TTL syntax is valid
- [ ] Memory management is correct
- [ ] Error handling implemented
- [ ] Crossfade logic preserved for smooth transitions
- [ ] Control port safety rules applied (null checks with defaults)
- [ ] Audio port validation present
- [ ] Latency compensation handles startup period correctly
- [ ] Processing called unconditionally, output routing controlled separately

## Code Organization

```
noise-repellent/
├── plugins/              # LV2 plugin implementations
├── lv2ttl/              # TTL file templates
├── src/                 # Shared utilities
├── subprojects/         # Dependencies (libspecbleach)
├── meson.build         # Build configuration
├── meson_options.txt   # Build options
```

## Performance Considerations
- Low latency processing (< 10ms typical)
- Efficient STFT implementation via libspecbleach
- Memory-aligned allocations
- Smooth parameter transitions via crossfade
