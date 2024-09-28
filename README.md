# noise-repellent-DEV-Efenstor-fork

### This is the fixed development branch of the original noise-repellent, probably not quite stable but tested and useable.

A suite of lv2 plugins for noise reduction that uses [libspecbleach](https://github.com/lucianodato/libspecbleach) C library.

## Features

* Adaptive noise reduction plugin for low latency voice denoise
* Manual noise capture based plugin for customizable noise reduction
* Adjustable Reduction and many other parameters to tweak the reduction
* Option to listen to the residual signal
* Soft bypass
* Noise profile saved with the session

## Install

To compile yourself and install this plug-in you will need the a C compiling toolchain, LV2 SDK, Meson build system, ninja compiler, git and libspecbleach library (if it doesn't find it it will download and compile it. In this case make sure to have libspecbleach dependencies installed).

Installation:

```bash
  git clone https://github.com/Efenstor/noise-repellent-DEV-Efenstor-fork.git --branch=DEVELOPMENT
  cd noise-repellent-DEV-Efenstor-fork
  mkdir build
  cd build
  meson build .. --buildtype=release
  meson compile -C build -v
  sudo meson install -C build
```

## Use Instuctions

Please refer to project's wiki <https://github.com/lucianodato/noise-repellent/wiki>
