#!/bin/bash

# Clean
rm -rf build lv2 nrepel.lv2 /tmp/lv2/nrepel.lv2 noise-repellent.zip || true

# Get lv2 ttls
git clone https://gitlab.com/lv2/lv2.git

# Build now
meson build --buildtype release --prefix /tmp --libdir /tmp
ninja -v -C build

# Run tests
ninja -C build test