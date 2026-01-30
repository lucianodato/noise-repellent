---
description: "Workflow for local development of libspecbleach integrated with noise-repellent"
---

# libspecbleach Integration Workflow

Use this workflow when making changes to `libspecbleach` that need to be tested within `noise-repellent`.

## 1. Setup Local Development

1. Create a symlink to your local `libspecbleach` repo:
   ```bash
   cd subprojects
   ln -sf ../../../libspecbleach libspecbleach
   ```
2. Configure with `forcefallback` to ensure the local version is used:
   ```bash
   meson setup build --wrap-mode=forcefallback
   ```
3. Verify the symlink is being used:
   ```bash
   meson compile -C build -v | grep libspecbleach
   ```

## 2. Iterative Development

1. Make changes to `libspecbleach`.
2. Format `libspecbleach` code:
   ```bash
   ninja -C libspecbleach/build format
   ```
3. Rebuild `noise-repellent`:
   ```bash
   meson compile -C build
   ```
4. Test the plugin functionality in a host or DAW.

## 3. Troubleshooting

- **API Changes**: If function signatures differ, update the calls in `noise-repellent`.
- **Linker Errors**: Check library versions and include paths.
- **Reconfiguration**: If things get messy, run:
  ```bash
  rm -rf build && meson setup build --wrap-mode=forcefallback
  ```
