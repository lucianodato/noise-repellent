# noise-repellent

A suite of lv2 plugins for noise reduction that uses [libspecbleach](https://github.com/lucianodato/libspecbleach) C library.

[![build](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml)

## Features

* Adaptive noise reduction plugin for low latency voice denoise
* Manual noise capture based plugin for customizable noise reduction
* Adjustable Reduction and many other parameters to tweak the reduction
* Option to listen to the residual signal
* Soft bypass
* Noise profile saved with the session

## Install

Binaries for most platforms are provided with Github release. Just extract the adequate zip file for your platform to your [lv2 plugins folder](https://lv2plug.in/pages/filesystem-hierarchy-standard.html)

If you wish to compile yourself and install this plug-in you will need the a C compiling toolchain, LV2 SDK, Meson build system, ninja compiler, git and libspecbleach library (if it doesn't find it it will download and compile it. In this case make sure to have libspecbleach dependencies installed). It is recommended for you to include some optimizations in environment CFLAGS if you are in x86_64 architecture like the suggested bellow but you can customize them for the architecture that you are in.

Installation:

```bash
  git clone https://github.com/lucianodato/noise-repellent.git
  cd noise-repellent
  CFLAGS="-ffast-math -msse -msse2 -mfpmath=sse" meson build --buildtype=release --prefix=/usr --libdir=lib (your-os-appropriate-location-fullpath)
  meson compile -C build -v
  sudo meson install -C build
```

Noise-repellent is on Arch community at <https://www.archlinux.org/packages/community/x86_64/noise-repellent/>.

Noise-repellent is also available in KXStudios repositories <https://kx.studio/Repositories:Plugins>

## Use Instuctions

Please refer to project's wiki <https://github.com/lucianodato/noise-repellent/wiki>
