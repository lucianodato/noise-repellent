#!/bin/bash

mkdir -p build
cd build
meson .. --buildtype release --strip
ninja -j2 -v
sudo ninja install


