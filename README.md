noise-repellent
-------
An lv2 plugin for broadband noise reduction.

Features
-------
* Spectral gating and spectral sustraction supression rule
* Adaptive and manual noise threshold estimation
* Time smoothing and envelopes for musical noise reduction
* Regulable noise floor
* Regulable offset of thresholds to perform oversustraction
* Wideband Post-gate for cleaning noise only zones
* Whitening of the noise floor
* Option to listen to the residual signal
* Soft bypass
* Noise print saved with the session

Limitations
-------
* The plugin will introduce latency so it's not appropriate to be used while recording
* It was developed to be used with Ardour in mind  however it is known to work with other host

Install
-------
To compile and install this plugin you will need the LV2 SDK, gnu-make, a c-compiler, git, pkg-config and fftw3 library.

```bash
  git clone https://github.com/lucianodato/noise-repellent.git
  make
  sudo make install
```
Usage
-----
Manual noise capture workflow:
1) First select a portion of noise in your track of at least one second and loop it. The longer the better.
2) Turn on Capture noise print for a bit (at least one loop).
3) Once that's done turn it off and tweak parameters as you like.
Adaptive noise capture workflow:
1) Turn on Adpative mode and keep it on
2) Tweak parameters (artifact control and smoothing won't work)


Control Ports explained
-----
* Amount of reduction: Determines how much the noise floor will be reduced.
* Noise Offset: Scales the noise print captured. Greater values will reduce more noise at the expense of removing low level detail of the signal.
* Release: Timing of spectral gate releases. Larger values will reduce artifacts but might blur nearby transients.
* Smoothing: Reduces the variance of the reduction between frames to remove musical noise. Greater values may introduce echoes in the signal and will weaken transients.
* Artifact control: Interpolate between spectral gating and wideband gating in low SNR (pauses) zones. Greater values might reduce artifacts in those zones.
* Whitening: Modifies the residual noise to be more like white noise. This takes into account that our ears do well discriminating sounds in white noise versus colored noise.
* Makeup Gain: Output gain if needed.
* Capture noise print: To manually take the noise print.
* Adaptive Noise: To change the noise profile dynamically in time. This enables the automatic estimation of noise thresholds. It needs a few seconds to learn it.
* Reset noise print: Removes the noise print previously captured.
* Noise listen: To hear only the residual noise.
* Residual emphasis: When whitening is applied this will reduce lower frequencies but keep the boosted highs.
