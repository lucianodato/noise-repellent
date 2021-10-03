#!/bin/bash

# Clean
rm -rf build nrepel.lv2 /tmp/lv2/nrepel.lv2 noise-repellent.zip || true

# Build now
meson build --buildtype release --prefix /tmp --libdir /tmp
ninja -v -C build
ninja -C build install

# Compress build in a zip file
mv /tmp/lv2/nrepel.lv2 ./nrepel.lv2
zip -9 -r noise-repellent.zip ./nrepel.lv2