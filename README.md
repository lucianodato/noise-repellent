# noise-repellent

An suite of lv2 plugins for noise reduction that uses libspecbleach.

[![build](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml)

## Demo

[![Demo 1](http://img.youtube.com/vi/iNVxCvgcnig/0.jpg)](http://www.youtube.com/watch?v=iNVxCvgcnig "")
[![Demo 2](http://img.youtube.com/vi/LeKyGoAmbFE/0.jpg)](https://www.youtube.com/watch?v=LeKyGoAmbFE "")

## Features

* Adaptive noise reduction plugin for low latency voice denoise
* Manual noise capture plugin for customizable noise reduction
* Adjustable Reduction and many other parameters to tweak the reduction
* Option to listen to the residual signal
* Soft bypass
* Noise profile saved with the session

## Limitations

* The plug-in will introduce latency so it's not appropriate to be used while recording (35 ms for 44.1 kHz)
* It was developed to be used with Ardour however it is known to work with other hosts

## Install

Binaries for most platforms are provided with Github release. Just extract the adequate zip file for your platform to your lv2 plugins folder (normally /usr/local/lib/lv2 or $HOME/.lv2)

If you wish to compile yourself and install this plug-in you will need the LV2 SDK, Meson build system, ninja compiler, git and fftw3 library.

Installation:

```bash
  git clone https://github.com/lucianodato/noise-repellent.git
  cd noise-repellent
  meson build --buildtype=release --prefix=/usr --libdir=lib (your-os-appropriate-location-fullpath)
  ninja -C build -v
  sudo ninja -C build install
```

Noise-repellent is on Arch community at <https://www.archlinux.org/packages/community/x86_64/noise-repellent/>.

Noise-repellent is also available in KXStudios repositories <https://kx.studio/Repositories:Plugins>

## Use Instuctions

Please refer to project's wiki <https://github.com/lucianodato/noise-repellent/wiki>
