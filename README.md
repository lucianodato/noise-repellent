# noise-repellent

An lv2 spectral noise reduction plugin.

For Debian based systems you will need this packages to compile: lv2-dev libfftw3-dev

Compiling instructions: make and sudo make install

This plugin is intended to be used with Ardour. It will introduce latency so it only can be used on tracks in Ardour and not in busses.
It's a work in progress so don't expect it to be like the best comercial options.

Instructions: Select a section of noise only in your track and loop it. Turn on noise capture for one second or two to learn the noise profile and then turn it off. You can now adjust the reduction. For Hum noise use Max with CMSR. For White noise use Spectral Subtraction with Average and smoothing as you like. Much higher reduction strenght will do the job with the remaining musical noise..

All the parameters are exposed right now, so it may be a little confusing at first. I will be changing defaults values as I experiment with it and make it easier to use.

What is under the hood:

CMSR: It's a modified Ephraim-Malah noise supression rule made by Canazza-Mian to prevent transient supression and to be less echoie.

Spectral Subtraction: It's the classic power Subtraction supression rule.

Both Wiener and Ephraim-Malah algorithms are for testing purposes.

CMSR is supposed to introduce less musical noise, but it will blur your transients a bit (much less than the Ephraim-Malah supression though).

Spectral Substraction will keep transients "safe" but it will introduce musical noise. You coud get somewhat the same results as CMSR by using more detection smoothing.

Experiment with both they are both good.

Frequency smoothing uses a Savinzky-Golay quadratic filter. Other smoothings are just exponetial smoothing done between frames.
