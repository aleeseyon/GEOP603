"""Tests of ffd2d against exact/analytic limits. Run: pytest -q notebooks/tests"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import ffd2d as f  # noqa: E402


def test_fd_coeffs_known_limits():
    # 2nd order: [-2, 1]; 4th order: [-30, 16, -1]/12 (Taylor exactness)
    assert np.allclose(f.fd_coeffs(2), [-2.0, 1.0])
    assert np.allclose(f.fd_coeffs(4) * 12, [-30.0, 16.0, -1.0])
    # exact on x^2 and x^4 (4th order): sum c_m m^2 * 2 = 2, sum c_m m^4 = 0
    c = f.fd_coeffs(4)
    m = np.arange(len(c))
    assert np.isclose(2 * np.sum(c[1:] * m[1:] ** 4), 0.0)


def test_ffd_coeffs_reduce_to_identity_when_v_equals_v0():
    a, bx, bz = f.ffd_coeffs(2000.0, 2000.0, 1e-3, 10.0, 10.0)
    assert np.isclose(a, 1.0) and bx == 0.0 and bz == 0.0


@pytest.mark.parametrize("r", [0.6, 1.3, 2.0])
def test_ffd_symbol_matches_exact_to_4th_order_at_small_k(r):
    # S_ffd - S_exact = O(k^6): ratio of errors at k and k/2 ~ 2^6 = 64
    v, dt, dx = 2000.0, 1e-3, 10.0
    th = 0.3
    errs = []
    for k in (0.02, 0.01):
        kx, kz = k * np.cos(th), k * np.sin(th)
        e = f.symbol("ffd", kx, kz, v, dt, dx, dx, v / r) - f.symbol("exact", kx, kz, v, dt, dx, dx)
        errs.append(abs(e))
    assert 40 < errs[0] / errs[1] < 90


def test_green_freq_vs_time_domain():
    t = np.arange(0, 0.6, 1e-3)
    a = f.green2d_trace(500.0, 2000.0, t, 15.0, 0.1)
    b = f.green2d_trace_time(500.0, 2000.0, t, 15.0, 0.1)
    assert np.abs(a - b).max() / np.abs(b).max() < 2e-3


def test_ffd_homogeneous_v0_equals_v_is_exact_for_plane_wave():
    # free (no source) evolution of a Fourier mode must be cos(w t) exactly
    n, dx, v, dt = 64, 10.0, 2000.0, 2e-3  # Courant 0.4
    x = np.arange(n) * dx
    kx = 2 * np.pi * 5 / (n * dx)
    p1 = np.tile(np.cos(kx * x), (n, 1))
    w = v * kx
    p0 = np.tile(np.cos(kx * x), (n, 1)) * np.cos(w * dt)  # p(t=-dt)
    V = np.full((n, n), v)
    L = f.make_operator("ffd", V, dt, dx, dx, v0=v)
    pm, p = p0, p1
    for _ in range(200):
        pm, p = p, 2 * p - pm + L(p)
    exact = np.tile(np.cos(kx * x), (n, 1)) * np.cos(w * 200 * dt)
    assert np.abs(p - exact).max() < 1e-9


def test_ffd_homogeneous_matches_analytic_green():
    n, dx, v, dt, nt, f0, t0 = 256, 10.0, 2000.0, 1e-3, 450, 15.0, 0.1
    t = np.arange(nt) * dt
    V = np.full((n, n), v)
    tr, _, _ = f.propagate("ffd", V, dx, dx, dt, nt, (128, 128),
                           f.ricker(t, f0, t0), [(128, 178)], v0=v)
    an = f.green2d_trace(500.0, v, t, f0, t0)
    assert np.linalg.norm(tr[0] - an) / np.linalg.norm(an) < 5e-3


def test_fd2_fd4_stability_limits_2d():
    # analytic von Neumann limits for 2nd-order time, square grid, 2D:
    # FD2: max|S|/C^2 = 8 -> C <= 1/sqrt(2); FD4: 32/3 -> C <= sqrt(3/8)
    assert abs(f.max_stable_courant("fd2", 2000.0, 10.0, 10.0, n=1500) - 1 / np.sqrt(2)) < 3e-3
    assert abs(f.max_stable_courant("fd4", 2000.0, 10.0, 10.0, n=1500) - np.sqrt(3 / 8)) < 3e-3


def test_ffd_v0_equals_v_unconditionally_stable():
    for C in (0.5, 1.0, 2.0, 5.0):
        assert f.is_stable("ffd", 2000.0, C * 10.0 / 2000.0, 10.0, 10.0, v0=2000.0)


def test_predicted_ffd_limit_matches_time_stepping():
    # homogeneous v = 2 v0: predicted limit C*; run random field just below/above
    v, v0, dx, n = 2000.0, 1000.0, 10.0, 64
    Cs = f.max_stable_courant("ffd", v, dx, dx, v0=v0, n=600, nk=101)
    rng = np.random.default_rng(0)
    p0 = rng.standard_normal((n, n))
    V = np.full((n, n), v)
    growth = {}
    for C in (0.95 * Cs, 1.05 * Cs):
        dt = C * dx / v
        L = f.make_operator("ffd", V, dt, dx, dx, v0=v0)
        pm, p = p0.copy(), p0.copy()
        for _ in range(400):
            pm, p = p, 2 * p - pm + L(p)
        growth[C] = np.abs(p).max() / np.abs(p0).max()
    lo, hi = sorted(growth)
    assert growth[lo] < 10 and growth[hi] > 1e3
