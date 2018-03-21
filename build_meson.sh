#!/bin/bash

rm -r builddir
mkdir builddir
cd builddir
meson .. --buildtype release
ninja -v

