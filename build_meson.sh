#!/bin/bash

rm -r build
mkdir build
cd build
meson .. --buildtype debug
ninja -v

