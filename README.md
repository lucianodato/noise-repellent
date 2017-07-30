noise-repellent
-------
An lv2 plugin for broadband noise reduction.

Features
-------
* Spectral gating and spectral sustraction supression rule
* Adaptive and manual noise thresholds estimation
* Regulable noise floor
* Regulable offset of thresholds to perform oversustraction
* Time smoothing and a postfilter for musical noise reduction
* Whitening of the noise floor to mask musical noise
* Option to listen to the residual signal
* Soft bypass
* Noise profile saved with the session

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
Manual noise learn workflow:
1) First select a portion of noise in your track of at least one second and loop it. The longer the better.
2) Turn on learn noise profile for a bit (at least one loop).
3) Once that's done turn it off and tweak parameters as you like.

Adaptive noise learn workflow:
1) Turn on Adpative mode and keep it on
2) Tweak parameters (SNR threshold won't work)


Control Ports explained
-----
* Amount of reduction: Determines how much the noise floor will be reduced.
* Thresholds Offset: Scales the noise profile learned. Greater values will reduce more noise at the expense of removing low level detail of the signal.
* Release: Timing of the reduction applied. Larger values will reduce artifacts but might blur nearby transients.
* PF threshold: Threshold for the low SNR level detector of the postfilter. Lower thresholds will not apply the postfilter higher thresholds will start to discriminate between low and high SNR signals and apply filter smoothing to the lower ones.
* Whitening: Modifies the residual noise to be more like white noise. This takes into account that our ears do well discriminating sounds in white noise versus colored noise. Higher values will brighten the residual noise and will mask high frequency musical noise.
* Learn noise profile: To manually learn the noise profile.
* Adaptive Noise: To change the noise profile dynamically in time. This enables the automatic estimation of noise thresholds.
* Reset noise profile: Removes the noise profile previously learned.
* Residual listen: To hear only the residual noise.

Advice for better reduction
-----
General noise reduction advice
* Try to reduce and not remove entirely the noise. It will sound better
* Gentler settings with multiple intances of the plugin will probably sound better than too much reduction with one intance. Of course for every intance the noise should be re-learned again.
* If the noise varies to much from one section to other apply different reduction for each part.
* Always remember to listen to the residual signal to make sure that you are not distorting the signal too much.
* Start with the reduction at 0 dB and start decreasing it until you start to hear artifacts then tune the paremeters to get rid of them without distorting the signal.

For adaptive Reduction:
* Adaptive mode should be used only with voice tracks because the algorithm for noise estimation is tuned for that use.
* It's recommended to play with thresholds offset because those are estimated continuosly and the algorithm used for that tends to overestimate them.
* Shorter release will preserve higher frequencies better, but make sure to not set it so low that musical noise start to creep in.
* Set PF threshold to -60 dB to deactivate posffiltering. Postfilering will not do much when noise thresholds are being estimated continuosly.
* If the track you are processing does not have a long section of noise before the wanted signal starts cut some noise from an inbetween section and extend the beggining a bit. This is to take into account the time that takes the algorithm to learn the noise. Alternatively you can learn the noise by using one section of the track and then turning off the adaptive mode so a fixed noise profile is used. This will not adapt in time but will give you something to work with.

For manual reduction:
* You can use adaptive mode to estimate noise thresholds if there is no section in the track that contains noise only. It should be used the same way you would use the manual learn.
* If noise floor change a bit over time it might be useful to use higher thresholds offset.
* Make sure that the section you select to learn the noise profile is noise only (without breaths or sustained notes or anything but noise)
* The longer the section you select to learn the noise profile the better the reduction will sound.
* It's better to start reducing musical noise by using a longer release until some of the original signal is starting to appear in the residual noise. Reduce the release a bit and then increase the postfilter threshold for even more musical noise reduction. Of course check the residual noise for any distortion.
