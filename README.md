[![build](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml/badge.svg)](https://github.com/lucianodato/noise-repellent/actions/workflows/build.yml)

noise-repellent
-------
An lv2 plug-in for broadband noise reduction.

Demo
-------
[![](http://img.youtube.com/vi/iNVxCvgcnig/0.jpg)](http://www.youtube.com/watch?v=iNVxCvgcnig "")
[![](http://img.youtube.com/vi/LeKyGoAmbFE/0.jpg)](https://www.youtube.com/watch?v=LeKyGoAmbFE "")

Features
-------
* Spectral gating and spectral subtraction suppression rule
* Adaptive and manual noise thresholds estimation
* Adjustable noise floor
* Adjustable offset of thresholds to perform over-subtraction
* Time smoothing and a masking estimation to reduce artifacts
* Basic onset detector to avoid transients suppression
* Whitening of the noise floor to mask artifacts and to recover higher frequencies
* Option to listen to the residual signal
* Soft bypass
* Noise profile saved with the session

Limitations
-------
* The plug-in will introduce latency so it's not appropriate to be used while recording (35 ms for 44.1 kHz)
* It was developed to be used with Ardour however it is known to work with other hosts

Install
-------
Binaries for most platforms are provided with releases but if you are an experienced user you can go ahead an compile it from source. Just extract the adequate zip file for your platform to your lv2 plugins folder (normally /usr/local/lib/lv2 or $HOME/.lv2)

To compile and install this plug-in you will need the LV2 SDK, Meson build system (use pip3 to install it), ninja compiler, git and fftw3 library (>= 3.3.5 is recommended to avoid threading issues).

Installation:
```bash
  git clone https://github.com/lucianodato/noise-repellent.git
  cd noise-repellent
  meson build --buildtype release --prefix (your-os-appropriate-location-fullpath)
  ninja -v -C build
  sudo ninja -C build install
```
Noise-repellent is on Arch community at https://www.archlinux.org/packages/community/x86_64/noise-repellent/.

Noise-repellent is also available in KXStudios repositories https://kx.studio/Repositories:Plugins

Usage Instuctions
-----
Please refer to project's wiki https://github.com/lucianodato/noise-repellent/wiki
