# Noise Repellent — UX & Architectural Invariants

This document defines the critical UX invariants and architectural contracts that must be preserved across future refactors, optimizations, and feature additions.

---

## 1. Profile as Single Source of Truth
* **Persistence Across States:** Once noise is learned (`learn_noise`), the noise profile is frozen and preserved across playback starts, stops, transport suspensions, and smoothing mode switches. Switching smoothing mode (Standard (Fast & Low CPU) ↔ Patch-Based (High Quality) ↔ Patch-Based + Refinement (Max Quality), internally 1D Spectral ↔ 2D NLM ↔ 2D NLM + DFTT) must never lose the profile or cause the visualizer to flatline to -120 dB.
* **Unified Engine:** A single libspecbleach stereo group of unified spectral denoisers backs all modes; the smoothing strategy is selected per-block through parameters (`smoothing_mode`). The library owns the seamless mode transition (allocation-free internal crossfade, shared profile state; instant flip within the NLM family), so profile persistence across mode switches is structural, not synchronized.

## 2. Real-Time Interactivity While Silent / Stopped
* **Live Aggressiveness Control:** Adjusting the `aggressiveness` slider while playback is stopped or silent must immediately recalculate and display the morphed profile curve.
* **Live Threshold Offsets (Broadband & Tonal):** Adjusting `noise_profile_offset` or `tonal_noise_profile_offset` (when unlinked via `link_threshold_offset = false`) must immediately scale the displayed profile in real time.
* **Tonal Mask Synthesis:** If the engine has not yet executed its DSP pipeline, a synthetic tonal mask is generated from detected tonal peaks so that unlinked tonal offsets remain visually responsive.

## 3. Tonal Peak & Line Invariants
* **Visibility Rule:** Tonal lines and detected peaks are rendered if and only if `(!isLinked || !isOffsetLinked) && hasNoiseProfile`.
* **State Preservation:** Peak frequencies and tonal lines remain visible across smoothing mode switches and silence.

## 4. Silence Segregation
* **Input / Output Drop:** When the input signal drops below silence threshold (`< 1e-5` / ~-100 dB) and learning is inactive, input and output spectra drop to -120 dB.
* **Profile Isolation:** The noise floor profile must *never* drop or clear on silence; it remains frozen and interactive.
* **Flush Before Sleep:** The engine sleeps on silence only after it persists past the full pipeline (reported latency plus gain-release tails), so decaying tails ring out naturally and sparse transients are never swallowed mid-response.

## 5. Smoothing Mode Switching (Library-Owned)
* **Instantaneous From Plugin Perspective:** Switching between 1D Spectral, 2D NLM and 2D NLM + DFTT is a plain parameter load. libspecbleach performs the transition internally (allocation-free crossfade; instant flip within the NLM family since history and latency are shared); the plugin keeps a single engine instance and no switch machinery.
* **Constant Latency Per Frame Size:** For a given STFT frame size, all smoothing modes share the same algorithmic latency (temporal is padded to the NLM look-ahead; DFTT adds none). Changing the frame size (Options menu: 23 / 32 / 46 / 64 / 93 ms) changes latency (roughly 2x the frame) and re-reports PDC via a suspended rebuild; it must never rebuild concurrently with `processBlock`.
* **No Host Suspension For Mode Switches:** Mode switching must never trigger `suspendProcessing`, latency re-reports, or profile migration in the plugin. Frame-size switches are the exception: they suspend, rebuild, and re-report. Low-latency mode is the second exception: it suspends, rebuilds (512-sample frame, causal 1D-only, zero look-ahead), re-reports (~10.7 ms at 48 kHz / ~11.6 ms at 44.1 kHz), and starts clean (profile dropped, Learn auto-stopped) like a frame-size switch. While active, the smoothing selector is locked to Standard and the frame-size menu is disabled. The 512-sample frame is too coarse for independent tonal/broadband control, so reduction and threshold-offset links are forced on (unlinking unavailable).

## 7. Frame-Size Switching
* **Clean Slate:** A frame-size switch discards the learned profile (resampling across resolutions works poorly) and the HUD falls back to `NO PROFILE (PASS-THROUGH)` until the user re-learns at native resolution. Session state restores are exempt — profiles saved in a session load into the fresh engine normally.
* **Learn Safety:** An in-progress Learn is auto-stopped before the rebuild; a half-rolled mean must never migrate across resolutions.
* **Large-Frame Character (expected):** 64/93 ms trades time resolution for frequency resolution — transients smear across the (frame-counted) NLM patch and gain envelopes step per hop, which reads as "robotic" on transient material. This is inherent DSP tradeoff, not a defect; the Options menu steers transient material to 23/32 ms. Frame-rate-independent time constants are tracked as future library work (libspecbleach#152).

## 6. Real-Time Safety & Performance Hierarchy
* **Strict RT-Safety:** The audio thread (`processBlock`) must never perform dynamic memory allocations, take blocking locks/mutexes, or execute file/console I/O. Runtime smoothing mode changes are allocation-free inside the library.
* **Bypass Alignment:** The engine keeps running while bypassed (internal or host) so both ends of the mixer's crossfade carry pipeline-delayed content — toggling never time-travels by a full latency. The true-idle engine sleep applies only when unbypassed.
* **CPU Hierarchy:**
  $$\text{Silent/Idle (unbypassed)} < \text{Bypass} \approx \text{Denoise 1D} < \text{Denoise 2D} \le \text{Learn}$$
