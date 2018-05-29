[![Build Status](https://travis-ci.org/lucianodato/noise-repellent.svg?branch=master)](https://travis-ci.org/lucianodato/noise-repellent)

noise-repellent
-------
An lv2 plug-in for broadband noise reduction.

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
To compile and install this plug-in you will need the LV2 SDK, Meson build system (use pip3 to install it), ninja compiler, git and fftw3 library (>= 3.3.5 is recommended to avoid threading issues).

Installation:
```bash
  git clone https://github.com/lucianodato/noise-repellent.git
  cd noise-repellent
  chmod +x install.sh && ./install.sh
```
There is now an AUR package at https://aur.archlinux.org/packages/noise-repellent-git for Arch Users (Kindly done by CrocoDuck).

Code Documentation
-----
Code is documented using doxygen. To read it be sure to install doxygen in your system and run the following command:

```bash
  doxygen -s doc/doxygen.conf
```
This will generate an html folder inside doc folder. Accessing index.html you can read the documentation.