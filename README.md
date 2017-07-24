# noise-repellent

An lv2 spectral broadband noise reduction plugin.

Dependencies: Install your distribution packages similar to pkg-config lv2-dev and libfftw3-dev in Debian based systems.

Compiling instructions: make
Installing instructions: sudo make install

Warning!!! This plugin is meant to be used with Ardour. It will introduce latency. Don't expect it to be like commercial plugins, although it will do the job in most cases if configured well enough.

Using it: First select a portion of noise in your track of at least one second and loop it. The longer the better. Turn on Capture noise print for a bit (at least one loop). Use the reduction slider to reduce the noise.

Sliders explained:

- Amount of reduction: Determines how much the noise floor will be reduced.
- Noise Offset: Scales the noise print captured. Greater values will reduce more noise at the expense of removing low level detail of the signal.
- Release: Timing of spectral gate releases. Larger values will reduce artifacts but might blur nearby transients.
- Smoothing: Reduces the variance of the reduction between frames to remove musical noise. Greater values may introduce echoes in the signal and will weaken transients.
- Artifact control: Interpolate between spectral gating and wideband gating in low SNR (pauses) zones. Greater values might reduce artifacts in those zones.
- Whitening: Modifies the residual noise to be more like white noise. This takes into account that our ears do well discriminating sounds in white noise versus colored noise.
- Makeup Gain: Output gain if needed.

Buttons explained:

- Capture noise print: To manually take the noise print.
- Adaptive Noise: To change the noise profile dynamically in time. This enables the automatic estimation of noise thresholds. It needs a few seconds to learn it. (Smoothing and Artifact control won't work when this is active)
- Reset noise print: Removes the noise print previously captured.
- Noise listen: To hear only the residual noise.
- HF residual emphasis: When whitening is applied this will reduce lower frequencies but keep the boosted highs.
