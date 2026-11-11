# GEOP 603 term project plan

Plan, work queue, run log and decisions for the GEOP 603 term project. This
project is a separate objective from any other research; nothing here feeds
into, or is driven by, other work.

## Ground rules
1. No facts from recall: every equation, constant and claim about prior work
   is backed by a citation verified in the session, with its DOI (digital
   object identifier) in `report/refs.bib`. Unverifiable items are marked
   `UNVERIFIED:` and listed in DECISIONS_NEEDED.
2. All numbers come from executed code; notebooks are executed end-to-end
   (`jupyter nbconvert --execute`) before commit.
3. Each new function gets a test against an exact or analytic limit first.
4. Define every acronym at first use in every document.
5. Surgical changes only; never delete existing code (comment out and flag).
6. Numerical QC at the end of every item: stability (CFL,
   Courant–Friedrichs–Lewy), dispersion, units, signs.
7. Public repository: no personal information.

## Project

Song & Fomel (2011), Fourier finite-difference (FFD) wave propagation, Geophysics
76(5), T123–T129. Deliverables due December: executed notebooks, report,
20-min talk, defense Q&A.
Note (2026-11-10): the author's FFD files (BP 2007 TTI vz/vx/eta/theta .bin,
grid.json, `ffd_tti_devito.py`, notebooks) live in a separate interactive chat
on TTI dispersion. That chat concluded (sources: Song & Fomel 2011, Song,
Fomel & Ying 2013 GJI 193(2) doi:10.1093/gji/ggt017, Koene et al. 2018 GJI
213) that FFD dispersion is medium- and path-dependent, so no single post-hoc
Fourier filter removes it; in-loop fixes are lowrank FFD (LFFD), a
higher-order FD part, or several reference media. Re-verify those citations
before reuse. Do not duplicate that work; C3 may build on it once files arrive.
- [ ] C0. If the author's notebooks are in GEOP603: re-execute, log any
      failures, do not change physics. If absent: skip to C2/C5 and note it.
      (`notebooks/C1_verification.ipynb` and `ffd2d.py` were written by this
      queue in C1; they are not the author's notebooks and do not unblock C0.)
- [x] C1. Verification suite: homogeneous-medium FFD vs analytic 2D Green's
      function; dispersion vs 2nd/4th-order FD (finite difference) at equal
      Δt; stability limit scan vs reference-velocity choice.
- [ ] C2. Report skeleton `report/report.md` (Abstract, Theory with numbered
      equations from the paper, Implementation, Verification, Replication vs
      published figures, Extension, Conclusions); verified references.
- [ ] C3. Extension: **topic to be chosen by the author** (blocked). An
      earlier idea linked the extension to other research; it was removed on
      2026-11-11 so this project stays independent. See DECISIONS_NEEDED.
- [ ] C4. Slide outline + speaker notes for a 20-min talk (≈ 14 slides).
- [ ] C5. Defense Q&A bank: 30 likely examiner questions with short,
      referenced answers (stability, dispersion, reference velocity, TTI
      (tilted transverse isotropy), cost, limitations, extension).

Final interactive step (with the author): build the report as a Doc and the
talk as a Slides deck from the repo content.

## Run log
(one line per run: date, item, result, checks run, open issues)
- 2026-11-11 | setup | Plan moved here from a shared sprint file so that this course project is tracked separately | no code changed | Open: D1
- 2026-11-11 | C1 | Done: `notebooks/ffd2d.py` (2D FFD, FD2, FD4, symbols, von Neumann scan, analytic 2D Green's function), `notebooks/C1_verification.ipynb` (executed), 11 pytest tests, 4 figures in `report/figures/c1_*.png`. Homogeneous v=2000 m/s, dx=10 m, 15 Hz Ricker, r=1000 m: rel. L2 trace error vs analytic at dt=1 ms (Courant 0.2) FFD 0.18 % (v0=v), 0.19 % (v0=1500), 0.20 % (v0=2500); FD4 3.2 %; FD2 65.6 %. FFD error scales as dt^2 (ratios 4.02, 3.98 per halving). Phase velocity at 5 points/wavelength, dt=1 ms: FD2 0.938, FD4 0.990, FFD 0.9998-1.0003. Stability: FD2 C_max 0.706 (1/sqrt2), FD4 0.612 (sqrt(3/8)); FFD unconditional at v=v0, C_max 1.26-1.54 for v<v0 but 0.61 -> 0.50 for v/v0 = 1.05 -> 3 (confirmed by time stepping in tests). Paper's Fig. 3 set-up (dx=5 m, dt=2 ms, v0=RMS=801.5 m/s, v_max/v0=1.79, Courant 0.574) is unstable in this implementation (growth visible by 0.31 s); stable at dt=1 ms | Checks: pytest 11/11, nbconvert --execute, CFL table, dispersion curves, units/sign of source | Open: D2, D3


## DECISIONS_NEEDED
- D1 (2026-11-11): choose the C3 extension topic (course-only).
- D2 (2026-11-11, from C1): with eq. 8 as implemented (b_n = dt^2 v^2 (v^2 - v0^2)/(12 dx_n^2 v0^2), a = v^2/v0^2 - 2 sum b_n, re-derived by 2nd-order Taylor match and tested), the paper's smooth-model run (Fig. 3d: dx = 5 m, dt = 2 ms, v0 = RMS) blows up; the frozen-coefficient symbol at v_max reaches -4.645 < -4, only for |k| >= 1.29 x axis Nyquist (corners of the 2D wavenumber box). Possible causes: (a) the published code applies a high-wavenumber taper (the paper mentions suppressing high wavenumbers to improve stability, but the text does not specify one for Fig. 3); (b) a different v0 or dt in the actual run; (c) a transcription error in eq. 8 (the PDF text extraction was garbled; printed eq. 8 not yet checked by eye). Author: please check printed eq. 8 and, if available, the Madagascar script for Fig. 3; decide whether to add a radial |k| <= pi/dx taper (not implemented, to avoid speculative features).
- D3 (2026-11-11, from C1): the paper's smooth model with a 70 Hz maximum frequency at v_min = 550 m/s and dx = 5 m gives lambda_min = 7.9 m < 2 dx, i.e. |k| = 1.27 x Nyquist: under-sampled at the slowest velocity. Confirm what 'maximum frequency 70 Hz' means (peak vs maximum) before replicating Fig. 3; C1 used a 25 Hz peak Ricker, whose amplitude spectrum at 70 Hz is 0.84 % of its peak.
