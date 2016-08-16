# noise-repellent

An lv2 spectral noise reduction plugin.

For debian based systems you will need this packages to compile: lv2-dev lv2-c++-tools libfftw3-dev

Compiling instructions: make and make install then

This plugin is intended to be used with Ardour. The plugin introduce latency so it only can be used on tracks in Ardour and not in busses.
It's a work in progress so don't expect it to be like the best comercial options.

To use it you first have to select a section of noise in your track and play it looped
with the capture noise option ON. Once it's captured turn it off and play your track
while adjusting the reduction slider. If you hear tinkerbells make shure to play with the Reduce More slider
and with smoothing parameters to further reduce them. Use as less as posible with both to keep the inportant bits from the original signal.
Noise Whitening will recover some of the high frequency loss by filtering the residual noise.

Recomendations for now (This will be changed later as the dsp improves):
For Hum removal, a good choise is to use Max Noise Spectrum Statistic with the CMSR algorithm (Time Smoothing 3)
For White Noise like noise use Spectral subtraction algorithm with Geometric Mean statistic and More reduction slider
as far as it is not loosing too much high frequencies
