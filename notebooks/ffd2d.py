"""Minimal 2D acoustic Fourier finite-difference (FFD) propagator and checks.

Reference: Song, X., and S. Fomel, 2011, Fourier finite-difference wave
propagation: Geophysics, 76(5), T123-T129, doi:10.1190/geo2010-0287.1.

Wave equation (paper eq. 1, plus a point source added here):
    p_tt = v(x)^2 lap(p) + s(t) delta(x - xs)
Two-step time march (paper eq. 5):
    p(t+dt) + p(t-dt) - 2 p(t) = L p(t)
FFD operator (paper eqs. 6-9, 2D):
    L p = [a + 2 b_x cos(kx dx) + 2 b_z cos(kz dz)] applied (as a 5-point
          stencil in x) to IFFT{ 2[cos(v0 |k| dt) - 1] FFT{p} }
with a, b_n from the 2nd-order Taylor match about k = 0 (paper eq. 8):
    b_n = dt^2 v^2 (v^2 - v0^2) / (12 dx_n^2 v0^2)
    a   = v^2/v0^2 - 2 sum_n b_n
Reference FD schemes: 2nd-order time, 2nd/4th-order centred space.

Units: SI (m, s, m/s). Grid axis 0 = z, axis 1 = x. Boundaries periodic
(FFT and np.roll); domains are sized so wrap-around arrives after the
traces of interest.
"""
import numpy as np
from scipy.special import hankel2


# ---------------------------------------------------------------- sources
def ricker(t, f0, t0):
    """Ricker wavelet (second derivative of a Gaussian, unit peak) at t."""
    arg = (np.pi * f0 * (t - t0)) ** 2
    return (1.0 - 2.0 * arg) * np.exp(-arg)


# ---------------------------------------------------------------- coefficients
def fd_coeffs(order):
    """Centred 2nd-derivative weights c_0..c_M (stencil c_0 + sum c_m(u_{+m}+u_{-m}))/h^2.

    Solved from the Taylor (polynomial-exactness) conditions, which on an
    equispaced centred stencil gives the same weights as Fornberg (1988).
    """
    M = order // 2
    # unknowns c_0..c_M ; conditions: sum_m w_m m^(2j) = (2j)! delta_{j1}, j=0..M
    A = np.zeros((M + 1, M + 1))
    rhs = np.zeros(M + 1)
    for j in range(M + 1):
        A[j, 0] = 1.0 if j == 0 else 0.0
        for m in range(1, M + 1):
            A[j, m] = 2.0 * m ** (2 * j)
    rhs[1] = 2.0
    return np.linalg.solve(A, rhs)


def ffd_coeffs(v, v0, dt, dx, dz):
    """FFD coefficients a, b_x, b_z (Song & Fomel 2011, eq. 8, 2D form)."""
    v = np.asarray(v, float)
    bx = dt ** 2 * v ** 2 * (v ** 2 - v0 ** 2) / (12.0 * dx ** 2 * v0 ** 2)
    bz = dt ** 2 * v ** 2 * (v ** 2 - v0 ** 2) / (12.0 * dz ** 2 * v0 ** 2)
    a = v ** 2 / v0 ** 2 - 2.0 * (bx + bz)
    return a, bx, bz


# ---------------------------------------------------------------- symbols
def symbol(method, kx, kz, v, dt, dx, dz, v0=None):
    """Fourier symbol S(k) of the operator L, so that 2[cos(w dt)-1] = S(k)
    in a homogeneous medium. Exact (analytic) value: 2[cos(v|k|dt)-1]."""
    if method == "ffd":
        a, bx, bz = ffd_coeffs(v, v0, dt, dx, dz)
        kk = np.sqrt(kx ** 2 + kz ** 2)
        return 2.0 * (np.cos(v0 * kk * dt) - 1.0) * (
            a + 2 * bx * np.cos(kx * dx) + 2 * bz * np.cos(kz * dz))
    if method in ("fd2", "fd4"):
        c = fd_coeffs(int(method[2]))
        dxx = c[0] + sum(2 * c[m] * np.cos(m * kx * dx) for m in range(1, len(c)))
        dzz = c[0] + sum(2 * c[m] * np.cos(m * kz * dz) for m in range(1, len(c)))
        return v ** 2 * dt ** 2 * (dxx / dx ** 2 + dzz / dz ** 2)
    if method == "exact":
        return 2.0 * (np.cos(v * np.sqrt(kx ** 2 + kz ** 2) * dt) - 1.0)
    raise ValueError(method)


def phase_velocity_ratio(method, kx, kz, v, dt, dx, dz, v0=None):
    """Numerical/true phase velocity, from cos(w dt) = 1 + S/2. NaN if unstable."""
    S = symbol(method, kx, kz, v, dt, dx, dz, v0)
    c = 1.0 + S / 2.0
    with np.errstate(invalid="ignore"):
        w = np.where((c >= -1) & (c <= 1), np.arccos(np.clip(c, -1, 1)) / dt, np.nan)
    return w / (v * np.sqrt(kx ** 2 + kz ** 2))


def is_stable(method, v, dt, dx, dz, v0=None, nk=201, tol=1e-12):
    """von Neumann test for a homogeneous medium: stable iff -4 <= S(k) <= 0
    for all k in the Nyquist box (roots of z^2 - (2+S) z + 1 on |z| = 1)."""
    kx = np.linspace(-np.pi / dx, np.pi / dx, nk)
    kz = np.linspace(-np.pi / dz, np.pi / dz, nk)
    KX, KZ = np.meshgrid(kx, kz)
    S = symbol(method, KX, KZ, v, dt, dx, dz, v0)
    return bool(np.all(S <= tol) and np.all(S >= -4.0 - tol))


def max_stable_courant(method, v, dx, dz, v0=None, cmax=3.0, n=3001, nk=201):
    """Largest Courant number C = v dt/dx below which every C on a fine scan
    is stable (first failure on the scan; returns cmax if none)."""
    cs = np.linspace(cmax / n, cmax, n)
    last = 0.0
    for C in cs:
        if not is_stable(method, v, C * dx / v, dx, dz, v0, nk=nk):
            return last
        last = C
    return cmax


# ---------------------------------------------------------------- propagation
def _wavenumbers(nz, nx, dz, dx):
    kz = 2 * np.pi * np.fft.fftfreq(nz, dz)
    kx = 2 * np.pi * np.fft.fftfreq(nx, dx)
    KX, KZ = np.meshgrid(kx, kz)
    return np.sqrt(KX ** 2 + KZ ** 2)


def make_operator(method, v, dt, dx, dz, v0=None):
    """Return L(p) for one time step on a periodic grid; v may be 2D."""
    nz, nx = v.shape
    if method == "ffd":
        K = _wavenumbers(nz, nx, dz, dx)
        phase = 2.0 * (np.cos(v0 * K * dt) - 1.0)
        a, bx, bz = ffd_coeffs(v, v0, dt, dx, dz)

        def L(p):
            q = np.real(np.fft.ifft2(phase * np.fft.fft2(p)))
            return (a * q + bx * (np.roll(q, 1, 1) + np.roll(q, -1, 1))
                    + bz * (np.roll(q, 1, 0) + np.roll(q, -1, 0)))
        return L
    if method in ("fd2", "fd4"):
        c = fd_coeffs(int(method[2]))
        f = v ** 2 * dt ** 2

        def L(p):
            lap = c[0] * p * (1 / dx ** 2 + 1 / dz ** 2)
            for m in range(1, len(c)):
                lap = lap + c[m] * ((np.roll(p, m, 1) + np.roll(p, -m, 1)) / dx ** 2
                                    + (np.roll(p, m, 0) + np.roll(p, -m, 0)) / dz ** 2)
            return f * lap
        return L
    raise ValueError(method)


def propagate(method, v, dx, dz, dt, nt, src, wavelet, recs, v0=None,
              snap_it=None, p0=None, p1=None):
    """March p^{n+1} = 2p^n - p^{n-1} + L p^n + dt^2 s^n delta_h(xs).

    src: (iz, ix); recs: list of (iz, ix); wavelet: array len nt (s(t_n)).
    delta_h = 1/(dx dz) at the source node. Returns traces (nrec, nt),
    the snapshot at step snap_it (or None), and the max |p| history.
    """
    L = make_operator(method, v, dt, dx, dz, v0)
    pm = np.zeros(v.shape) if p0 is None else p0.copy()
    p = np.zeros(v.shape) if p1 is None else p1.copy()
    tr = np.zeros((len(recs), nt))
    amax = np.zeros(nt)
    snap = None
    iz, ix = src
    for it in range(nt):
        # record p at t_it (p holds time level it)
        for r, (rz, rx) in enumerate(recs):
            tr[r, it] = p[rz, rx]
        amax[it] = np.abs(p).max()
        if snap_it is not None and it == snap_it:
            snap = p.copy()
        pn = 2 * p - pm + L(p)
        pn[iz, ix] += dt ** 2 * wavelet[it] / (dx * dz)
        pm, p = p, pn
    return tr, snap, amax


# ---------------------------------------------------------------- analytic
def green2d_trace(r, v, t, f0, t0, pad=16):
    """Analytic p(r, t) for p_tt - v^2 lap p = s(t) delta(x), s = Ricker,
    homogeneous 2D, via frequency domain (time dependence e^{+i w t},
    numpy FFT convention):  P(w) = S(w) * (-i/4) H0^(2)(w r / v) / v^2.

    Normalisation follows from the 2D retarded Green's function
    G = H(vt - r) / (2 pi v sqrt(v^2 t^2 - r^2)) (Lousto & Nakano 2008 for
    v = 1, rescaled); test_green_freq_vs_time checks the two forms agree.
    The source is delayed by t0 and the record padded `pad`x to limit
    wrap-around of the 1/t tail. The discrete source term dt^2 s^n enters
    the second difference centred on t_n, so no extra delay is applied.
    """
    dt = t[1] - t[0]
    n = len(t) * pad
    tt = np.arange(n) * dt
    s = ricker(tt, f0, t0)
    S = np.fft.rfft(s)
    w = 2 * np.pi * np.fft.rfftfreq(n, dt)
    G = np.zeros_like(S)
    G[1:] = (-0.25j) * hankel2(0, w[1:] * r / v) / v ** 2
    p = np.fft.irfft(S * G, n)
    return p[: len(t)]


def green2d_trace_time(r, v, t, f0, t0):
    """Same as green2d_trace but by direct time-domain convolution with the
    retarded Green's function, using tau = r/v + u^2 to remove the 1/sqrt
    singularity (independent normalisation check; slow, for tests)."""
    from scipy.integrate import quad
    out = np.zeros(len(t))
    for i, ti in enumerate(t):
        T = ti - r / v
        if T <= 0:
            continue
        f = lambda u: ricker(ti - (r / v + u * u), f0, t0) / (
            np.pi * v * np.sqrt(v * (v * (r / v + u * u) + r)))
        out[i] = quad(f, 0.0, np.sqrt(T), limit=400)[0]
    return out
