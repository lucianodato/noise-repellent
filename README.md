# noise-repellent

An lv2 spectral noise reduction plugin.

For debian based systems you will need this packages to compile: lv2-dev lv2-c++-tools libfftw3-dev

Compiling instructions: make and make install then

This plugin is intended to be used with Ardour. The plugin introduce latency so it only can be used on tracks in Ardour and not in busses.
It's a work in progress so don't expect it to be like the best comercial options.

To use it you first have to select a section of noise in your track and play it looped
with the capture noise option ON. Once it's captured turn it off and play your track
while adjusting the reduction slider. If you hear tinkerbells make shure to play with the Over-sustraction slider
to further reduce them. Do not go overboard with it because it can distort your audio.
Noise Whitening will recover some of the high frequency loss by filtering the residual noise.
The rest of parameters are pretty safe with defaults.
