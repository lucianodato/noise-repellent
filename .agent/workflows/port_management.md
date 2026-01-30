---
description: "TTL file and port index synchronization checker and manager"
---

# Port Management Workflow

This workflow ensures that port indices in C code and TTL files remain synchronized when adding or removing ports.

## Adding a Port

1. Add entry to `PortIndex` enum in the plugin C file.
2. Add port definition to corresponding `.ttl.in` file in `lv2ttl/`.
3. Shift all subsequent port indices if inserting in the middle.
4. Rebuild to verify:
   ```bash
   meson compile -C build
   ```

## Removing a Port

1. Remove the entry from `PortIndex` enum in the C file.
2. Remove the port definition from the `.ttl.in` file.
3. Shift subsequent indices down by 1 in both files.
4. Update `connect_port()` switch cases in the C file.
5. Rebuild and test.

## Validation

Run these commands to verify synchronization:

```bash
# Check generated TTL files for indices
grep "lv2:index" build/*.ttl

# Verify C enum values
grep "NOISEREPELLENT_" plugins/*.c
```

---

## Example: Removing Enable Control

If removing a control like `NOISEREPELLENT_ENABLE`:

**C Code Changes:**
- Remove `NOISEREPELLENT_ENABLE` from `PortIndex` enum.
- Shift subsequent indices (e.g., LATENCY becomes 7, INPUT_1 becomes 8).
- Remove the `enable` field from the plugin struct.
- Remove the `enable` case from `connect_port()`.
- Update processing logic (e.g., change `signal_crossfade_run()` to always pass `true`).

**TTL Changes:**
- Remove the `enable` port definition (index 7).
- Shift subsequent ports: latency (7), input (8), output (9).
