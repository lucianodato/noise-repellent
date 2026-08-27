# Noise Repellent — UX & Architectural Invariants

This document defines the critical UX invariants and architectural contracts that must be preserved across future refactors, optimizations, and feature additions.

---

## 1. Profile as Single Source of Truth
* **Persistence Across States:** Once noise is learned (`learn_noise`), the noise profile is frozen and preserved across playback starts, stops, and transport suspensions.
* **Dual-Engine Coherence (1D ↔ 2D):** Profiles are synchronized between the 1D (Spectral) and 2D (NLM) engine groups. Switching algorithm mode while transport is stopped (or during playback) must never lose the profile or cause the visualizer to flatline to -120 dB.
* **Silent & First-Switch Fallback:** When playback is stopped or silent, the visualizer queries the active engine group; if the newly selected engine has not yet processed audio, it falls back to the previous group's profile and synthesizes necessary spectral state.

## 2. Real-Time Interactivity While Silent / Stopped
* **Live Aggressiveness Control:** Adjusting the `aggressiveness` slider while playback is stopped or silent must immediately recalculate and display the morphed profile curve.
* **Live Threshold Offsets (Broadband & Tonal):** Adjusting `noise_profile_offset` or `tonal_noise_profile_offset` (when unlinked via `link_threshold_offset = false`) must immediately scale the displayed profile in real time.
* **Tonal Mask Synthesis:** If an engine has not yet executed its DSP pipeline (e.g. immediately after a silent engine switch), a synthetic tonal mask is generated from detected tonal peaks so that unlinked tonal offsets remain visually responsive.

## 3. Tonal Peak & Line Invariants
* **Visibility Rule:** Tonal lines and detected peaks are rendered if and only if `(!isLinked || !isOffsetLinked) && hasNoiseProfile`.
* **State Preservation:** Peak frequencies and tonal lines remain visible across engine switches and silence.

## 4. Silence Segregation
* **Input / Output Drop:** When the input signal drops below silence threshold (`< 1e-5` / ~-100 dB) and learning is inactive, input and output spectra drop to -120 dB.
* **Profile Isolation:** The noise floor profile must *never* drop or clear on silence; it remains frozen and interactive.

## 5. Engine Switch Progress & Gapless Transitions
* **Gapless Phase Machine:** Engine switching proceeds through `Warming` (700 ms) → `XFade` (30 ms) → `Steady`.
* **Host Suspension Tolerance:** `getEngineSwitchProgress()` (0.0 → 1.0) must advance both during active `processBlock()` calls and via the processor's internal timer (`timerCallback` at 60 Hz) when the DAW transport is stopped or suspended. A silent switch completes in ~730 ms.

## 6. Real-Time Safety & Performance Hierarchy
* **Strict RT-Safety:** The audio thread (`processBlock`) must never perform dynamic memory allocations, take blocking locks/mutexes, or execute file/console I/O.
* **CPU Hierarchy:**
  $$\text{Bypass} < \text{Silent/Idle} < \text{Learn} < \text{Denoise 1D} < \text{Denoise 2D}$$
