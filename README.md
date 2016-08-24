# noise-repellent

An lv2 spectral noise reduction plugin.

Dependencies: Install your distribution packages that contain lv2-headers fftw3 libraries (lv2-dev and libfftw3-dev in Debian based systems case)

Compiling instructions: make and sudo make install

Warning!!! This plugin ment to be used with Ardour. It will introduce latency.

Using it: First select a portion of noise in your track of at least one second and loop it. Turn on Noise print capture for a bit (at least one loop). Use the reduction slider to reduce the noise.

How to improve the reduction quality: The strenght control will reduce more noise at the expense of removing low level detail in the signal. It wil sound less bright. The smoothing control will control the amount of musical noise that is present after the reduction, if it is pushed to hard will blur transients and sound echoie. A whitening option for post processing is provided where the residual spectrum will be modified to sound more like white noise, taking into consideration that our ears do well in discriminating sounds in white noise versus colored noise. Furthermore it will recover some of the high frequency detail and hide the low frequency noise more. Whitening is dependant of the reduction if you reduce a lot it won't have any effect since the noise floor is too low. A noise listen control is provided to tune what is reduced by listenig to the residual signal. Other thing that is important is that you select a representative noise to capture, since reduction is based on that. The longest the better. Make sure that is noise only too!
