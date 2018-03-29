#!/bin/bash

#remove previous build
rm -rf build || true
#remove previous install
sudo rm -rf /usr/local/lib/lv2/nrepel.lv2 || true

#make a new build directory
mkdir -p build

#build the plugin in the new directory
cd build
meson .. --buildtype release --strip
ninja -j2 -v

#install the plugin in the system
sudo ninja install


