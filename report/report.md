# Fourier finite-difference wave propagation: a replication of Song & Fomel (2011)

GEOP 603 term report — skeleton (item C2). Sections marked **TODO** are filled
by later queue items. Every number below comes from executed code in this
repository (`notebooks/`); every literature claim cites `report/refs.bib`.

Acronyms (defined at first use): FFD = Fourier finite difference; FD =
finite difference (FD2/FD4: 2nd/4th-order-in-space centred FD, 2nd order in
time); FFT/IFFT = (inverse) fast Fourier transform; RMS = root mean square;
CFL = Courant–Friedrichs–Lewy; RTM = reverse-time migration; TTI = tilted
transverse isotropy; LFFD = lowrank FFD; L2 = Euclidean (least-squares) norm.

## Abstract

**TODO** (after C3/replication). Draft scope: we re-derive and implement the
2D acoustic FFD time-stepping operator of Song & Fomel (2011)
[@SongFomel2011], verify it against the analytic 2D Green's function, compare
its dispersion and stability with FD2/FD4, and replicate the paper's
smooth-model experiment.

## 1. Theory

Equations are numbered as in Song & Fomel (2011) [@SongFomel2011]. Source of
the transcription: the author-hosted preprint PDF (ahay.org,
`RSF/book/tccs/ffd/paper.pdf`), read by text extraction on 2026-11-13.
**Caveat:** the printed journal version has not yet been compared by eye (see
DECISIONS_NEEDED D2(c) in `PLAN.md`).

Acoustic wave equation, pressure p, velocity v(x):

(1) ∂²p/∂t² = v(x)² ∇²p

Spatial Fourier transform, with p̂(k, t) = ∫ p(x, t) e^{ik·x} dx (eq. 3), in a
homogeneous medium:

(2) d²p̂/dt² = −v²|k|² p̂

(4) p̂(k, t + Δt) = e^{±i|k|vΔt} p̂(k, t)

Adding the two signs gives the two-step pseudo-analytical march

(5) p(x, t+Δt) + p(x, t−Δt) − 2p(x, t) = 2∫ p̂(k, t)[cos(|k|vΔt) − 1] e^{−ik·x} dk

which is exact for constant v. For v = v(x), factor out a constant reference
velocity v₀:

(6) 2[cos(v(x)|k|Δt) − 1] = 2[cos(v₀|k|Δt) − 1] · [cos(v(x)|k|Δt) − 1] / [cos(v₀|k|Δt) − 1]

and approximate the ratio by a short cosine (FD) series in each axis n with
grid spacing Δxₙ:

(7) [cos(v|k|Δt) − 1] / [cos(v₀|k|Δt) − 1] ≈ a + 2 Σₙ bₙ cos(kₙΔxₙ)

(8) a = (v²/v₀²)[1 + Δt²(v₀² − v²) Σₙ(1/Δxₙ²) / 6],  bₙ = Δt² v² (v² − v₀²) / (12 Δxₙ² v₀²)

(In 3D the paper writes Σₙ 1/Δxₙ² as (Δx₁²Δx₂² + Δx₂²Δx₃² + Δx₃²Δx₁²)/(Δx₁²Δx₂²Δx₃²).)
Equation 8 is equivalent to a = v²/v₀² − 2Σₙbₙ, i.e. a + 2Σbₙ = v²/v₀²
(the k → 0 match of eq. 7); this equivalence is tested in
`test_ffd_coeff_a_matches_paper_eq8_printed_form`. In practice the FFD step
is: FFT → multiply by 2[cos(v₀|k|Δt) − 1] → IFFT → apply the 3-point-per-axis
stencil with coefficients a, bₙ (eq. 9 in 3D; eq. 12 is the 4th-order
version with extra coefficients c, dₙ). Eqs. 10–11 give the TTI phase
function and rotated wavenumbers and are outside the acoustic scope of this
report unless the extension (C3) uses them.

FFD descends from Fourier finite-difference migration in the one-way
(depth-extrapolation) setting [@RistowRuhl1994]; lowrank variants (LFD,
LFFD) were developed later by the same group [@SongFomelYing2013].

## 2. Implementation

`notebooks/ffd2d.py` (2D, periodic boundaries, SI units, axis 0 = z):

- `ffd_coeffs` — eq. 8 in 2D; `make_operator('ffd', …)` — eqs. 5–7, one
  FFT/IFFT pair per step plus a 5-point stencil.
- FD2/FD4 reference schemes with Taylor-derived centred weights, identical
  to Fornberg's weights on an equispaced stencil [@Fornberg1988].
- Point source: dt² s(tₙ) δ_h added at the source node, δ_h = 1/(ΔxΔz);
  Ricker wavelet.
- `symbol`, `is_stable`, `max_stable_courant` — homogeneous-medium Fourier
  symbol and von Neumann test (stable iff −4 ≤ S(k) ≤ 0 over the Nyquist
  box) [@Lines1999].

**TODO:** absorbing boundaries are not implemented (not needed for C1).

## 3. Verification (item C1; `notebooks/C1_verification.ipynb`)

Homogeneous medium, v = 2000 m/s, Δx = Δz = 10 m, 15 Hz Ricker, receiver at
r = 1000 m. Reference: analytic 2D Green's function [@LoustoNakano2008],
computed in the frequency domain and cross-checked by direct time-domain
convolution (test `test_green_freq_vs_time_domain`).

Table 1. Relative L2 trace error (%) vs analytic, receiver at 0° (and 45°).

| Δt (ms) | Courant | FFD v₀=2000 | FFD v₀=1500 | FFD v₀=2500 | FD4 | FD2 |
|---|---|---|---|---|---|---|
| 1.0 | 0.20 | 0.182 (0.182) | 0.191 (0.181) | 0.204 (0.192) | 3.242 (1.813) | 65.644 (36.313) |
| 2.0 | 0.40 | 0.731 (0.731) | 0.768 (0.730) | 0.834 (0.778) | 7.205 (10.825) | 59.459 (27.598) |

- Temporal convergence (FFD, v₀ = v): error 0.7315 %, 0.1818 %, 0.0457 % at
  Δt = 2, 1, 0.5 ms; ratios 4.02, 3.98 → second order in Δt, as expected for
  the two-step march with a source sampled at tₙ (with v₀ = v the spatial
  operator is exact, so the residual is source discretisation; **TODO:**
  confirm this attribution with a band-limited initial-value test).
- Phase velocity ratio at 5 points per wavelength, Δt = 1 ms: FD2 0.938,
  FD4 0.990, FFD 0.9998–1.0003 for v₀ ∈ {1000, 1500, 2500} m/s.
- Stability (Courant number C = vΔt/Δx): FD2 C_max = 0.706 (theory 1/√2),
  FD4 0.612 (theory √(3/8)). FFD is unconditionally stable at v = v₀; for
  v < v₀, C_max = 1.26–1.54; for v > v₀, C_max falls from 0.610 (v/v₀ = 1.05)
  to 0.495 (v/v₀ = 3). Figures: `figures/c1_traces_vs_analytic.png`,
  `figures/c1_dispersion.png`, `figures/c1_stability_vs_v0.png`.

## 4. Replication of published figures

Smooth model (paper Fig. 3): v(x, z) = 550 + 1.5×10⁻⁴(x − 800)² +
10⁻⁴(z − 500)² m/s, 0 ≤ x, z ≤ 2560 m, Δx = 5 m, Δt = 2 ms, single reference
velocity v₀ = RMS velocity, Ricker source of "maximum frequency 70 Hz" at
the model centre [@SongFomel2011].

Status from C1: on the discrete grid (v range 550.0–1434.3 m/s; the paper
quotes 550–1439 m/s, the continuous maximum at x = z = 2560 m), v₀ = 801.5 m/s
and v_max/v₀ = 1.79. At Δt = 2 ms the frozen-coefficient symbol reaches
−4.645 < −4 for |k| ≥ 1.289 × axis Nyquist and the run grows without bound
from 0.308 s; at Δt = 1 ms it is stable. Open questions D2, D3 (`PLAN.md`).
Figure: `figures/c1_smooth_model_stability.png`.

**TODO:** BP model (paper Figs. 4–6) and TTI examples (Figs. 7–9) — data
availability not yet established.

## 5. Extension

**TODO** — topic to be chosen by the author (D1).

## 6. Conclusions

**TODO.**

## References

Keys refer to `report/refs.bib` (all DOIs checked against Crossref on
2026-11-13): SongFomel2011, RistowRuhl1994, SongFomelYing2013, Fornberg1988,
Lines1999, LoustoNakano2008.
