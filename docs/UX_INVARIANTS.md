# Noise Repellent — UX & Architectural Invariants

This document defines the critical UX invariants and architectural contracts that must be preserved across future refactors, optimizations, and feature additions.

---

## 1. Profile as Single Source of Truth
* **Persistence Across States:** Once noise is learned (`learn_noise`), the noise profile is frozen and preserved across playback starts, stops, transport suspensions, and smoothing mode switches. Switching smoothing mode (Standard (Fast & Low CPU) ↔ Patch-Based (High Quality), internally 1D Spectral ↔ 2D NLM) must never lose the profile or cause the visualizer to flatline to -120 dB.
* **Unified Engine:** A single libspecbleach stereo group of unified spectral denoisers backs both modes; the smoothing strategy is selected per-block through parameters (`smoothing_mode`). The library owns the seamless mode transition (allocation-free internal crossfade, shared profile state), so profile persistence across mode switches is structural, not synchronized.

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

## 5. Smoothing Mode Switching (Library-Owned)
* **Instantaneous From Plugin Perspective:** Switching between 1D Spectral and 2D NLM is a plain parameter load. libspecbleach performs the transition internally (allocation-free crossfade); the plugin keeps a single engine instance and no switch machinery.
* **Constant Latency:** Both smoothing modes share the same algorithmic latency (temporal is padded to the NLM look-ahead), so `getLatencySamples()` and host PDC never change on a mode switch.
* **No Host Suspension:** Mode switching must never trigger `suspendProcessing`, latency re-reports, or profile migration in the plugin.

## 6. Real-Time Safety & Performance Hierarchy
* **Strict RT-Safety:** The audio thread (`processBlock`) must never perform dynamic memory allocations, take blocking locks/mutexes, or execute file/console I/O. Runtime smoothing mode changes are allocation-free inside the library.
* **CPU Hierarchy:**
  $$\text{Bypass} < \text{Silent/Idle} < \text{Learn} < \text{Denoise 1D} < \text{Denoise 2D}$$
