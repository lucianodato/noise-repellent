# noise-repellent

An lv2 spectral noise reduction plugin.

Dependencies: Install your distribution packages similar to pkg-config lv2-dev and libfftw3-dev in Debian based systems.

Installing instructions: make and sudo make install

Warning!!! This plugin is meant to be used with Ardour. It will introduce latency. Don't expect it to be like commercial options, although it will do the job in most cases if configured well enough.

Using it: First select a portion of noise in your track of at least one second and loop it. Turn on Noise print capture for a bit (at least one loop). Use the reduction slider to reduce the noise.

How to improve the reduction quality:

- The strength control will reduce more noise at the expense of removing low level detail in the signal. It will sound less bright.

- The smoothing control will control the amount of musical noise that is present after the reduction, if it is pushed too hard it will blur transients and sound reverberant.

- NEW! Masking will use masking thresholds in order to estimate the strength adaptively and achieve much more noise reduction without distorting the signal too much and reducing musical noise too.

- A whitening option for post processing is provided where the residual spectrum will be modified to sound more like white noise, taking into consideration that our ears do well in discriminating sounds in white noise versus colored noise. Furthermore it will recover some of the high frequency detail and hide the low frequency noise more. Whitening is dependant of the reduction. If you reduce a lot, it won't have any effect since the noise floor is too low.

- A noise listen control is provided to tune what is reduced by listenig to the residual signal.

- Other thing that is important is that you select a representative noise region to capture, since reduction is based on that. The longest the better. Make sure that is noise only too!

- The last control available is frequency smoothing, it will reduce musical noise at the expense of some low frequency.

- An auto learn noise feature is provided. It tries to track the noise spectrum automaticaly at the same time signal is processed. It will work best for gentle noise reduction of speech signals. Use this feature where noise is changing over time and rapidly. It will never sound as good as if you do good noise profiling yourself.

General Rules:

- With hum noise, use strength and smoothing only.

- With white noise you can use frequency smoothing too. If you loose too much detail from the original signal, use whitening and not too much reduction. Multiple pass will reduce more noise in pauses and less during signal presence. If you use multiple passes, do gentle parameters each time and re learn noise every time.
