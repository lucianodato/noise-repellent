#!/bin/bash

rm -rf build || true
sudo rm -rf /usr/local/lib/lv2/nrepel.lv2 || true
mkdir -p build
cd build
meson .. --buildtype release --strip
ninja -j2 -v
sudo ninja install


