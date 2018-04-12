#!/bin/bash

#remove previous build
rm -rf build || true
#remove previous install
sudo rm -rf /usr/local/lib/lv2/nrepellent.lv2 || true

#build the plugin in the new directory
meson build --buildtype release
cd build
ninja -v

#install the plugin in the system
sudo ninja install


