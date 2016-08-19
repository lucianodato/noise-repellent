# noise-repellent

An lv2 spectral noise reduction plugin.

For Debian based systems you will need this packages to compile: lv2-dev lv2-c++-tools libfftw3-dev

Compiling instructions: make and make install

This plugin is intended to be used with Ardour. It will introduce latency so it only can be used on tracks in Ardour and not in busses.
It's a work in progress so don't expect it to be like the best comercial options.

Instructions: Select a section of noise only in your track and loop it. Turn on noise capture for one second or two to learn the noise profile and then turn it off. You can now adjust the reduction. For Hum noise use Max with CMSR. For White noise use Spectral Subtraction with Average and smoothing as you like and much higher reduction strenght.
