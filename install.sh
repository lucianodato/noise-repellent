#!/usr/bin/env bash

#default installation directories
INSTALL_DIR_LINUX="/usr/local/lib/lv2"
INSTALL_DIR_MAC="/Library/Audio/Plug-Ins/LV2"

# Detect the platform (similar to $OSTYPE)
OS="`uname`"
case $OS in
  'Linux') OS='Linux' && echo "You are on a Linux system. Building for Linux";;
  'Darwin') OS='Mac' && echo "You are on a Mac system. Building for MacOS";;
  *) ;;
esac

#remove previous builds
rm -rf build || true

#build the plugin in the new directory
if [ $OS = "Linux" ]; then
    meson build --buildtype release --prefix $INSTALL_DIR_LINUX
elif [ $OS = "Mac" ]; then
    meson build --buildtype release --prefix $INSTALL_DIR_MAC
fi

cd build
ninja -v

#install the plugin in the system
sudo ninja install