#!/bin/bash

#remove previous build
rm -rf build || true
sudo rm -rf /Library/Audio/Plug-Ins/LV2/nrepel.lv2 || true


#build the plugin in the new directory
meson build --buildtype release
cd build
ninja -v

#install the plugin in the system
sudo ninja install


