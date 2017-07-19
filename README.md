# noise-repellent

An lv2 spectral noise reduction plugin.

Dependencies: Install your distribution packages similar to pkg-config lv2-dev and libfftw3-dev in Debian based systems.

Installing instructions: make and sudo make install

Warning!!! This plugin is meant to be used with Ardour. It will introduce latency. Don't expect it to be like commercial options, although it will do the job in most cases if configured well enough.

Using it: First select a portion of noise in your track of at least one second and loop it. Turn on Noise print capture for a bit (at least one loop). Use the reduction slider to reduce the noise.

Sliders explained:

- Amount of reduction: Determines how much the noise floor will be reduced.
- Noise Offset: Scales the noise print captured. Greater values will reduce more noise at the expense of removing low level detail of the signal.
- Smoothing: Reduces the variance between frames to remove musical noise. Greater values may introduce echoes in the signal or blur transients.
- Artifact control: Interpolate between spectral gating and wideband gating in low SNR zones. Greater values might reduce artifacts in those zones.
- Release: Timing of spectral gate releases. Larger values will reduce artifacts but might blur nearby transients
- Whitening: Modifies the residual noise to be more like white noise. This take into account that our ears do well discriminating sounds in white noise versus colored noise.
- Makeup Gain: Output gain if needed.

Buttons explained:

- Capture noise print: To manually take the noise print.
- Adaptive Noise: To change the noise profile dynamically in time. This enables the automatic estimation of noise thresholds.
- Reset noise print: Removes the noise print previously captured.
- Noise listen: To hear only the residual noise.
- HF residual emphasis: Applies a window to the residue for HF boost. This might recover HF details.
- Transient preservation: Compute onsets to avoid distorting transients. Useful for percusive instruments.
