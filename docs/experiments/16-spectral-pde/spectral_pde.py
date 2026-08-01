#!/usr/bin/env python3
"""Experiment 16 -- Spectral PDE: an FFT inside a time loop (NumPy column).

Runs the same algorithm as ``spectral_pde.m``, in the same order, with the
same sizes.  See ``README.md`` for the measurements and the analysis.

    python3 spectral_pde.py

WHAT IT MEASURES.  The suite's existing spectral rows time ONE large
transform, so what they measure is FFTW.  A pseudo-spectral solver has the
opposite profile: thousands of MEDIUM transforms with complex elementwise
algebra between them, which measures per-call overhead and complex-array
arithmetic instead.

THE TRANSFORM CONVENTION IS LOAD-BEARING.  Wolfram's ``Fourier`` is
``(1/sqrt(n)) sum x e^{+2 pi i ...}`` -- a POSITIVE exponent -- which is
``np.fft.ifft(x, norm="ortho")``, and ``InverseFourier`` is
``np.fft.fft(x, norm="ortho")``.  The helpers below are named ``_F`` and
``_IF`` after the Wolfram functions they stand for, precisely so the mapping
is visible at every use.  Getting it backwards gives a plausible-looking wrong
answer -- the solution still evolves and still looks like turbulence -- and
only the cross-system check catches it.

THE CHECKS ARE SHORT RUNS ON PURPOSE.  Both equations are chaotic, so the
state at the final time is a property of the arithmetic rather than of the
equation.  50 KS steps and 20 NS steps are still inside the deterministic
regime.
"""

import time

import numpy as np

# ---- shared reporting helpers (identical in every experiment file) --------


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs."""
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))


def nest(f, x, m):
    for _ in range(m):
        x = f(x)
    return x


# Wolfram's Fourier / InverseFourier, exactly.  See the module docstring.
def _F(a):
    return np.fft.ifft(a, norm="ortho")


def _IF(a):
    return np.fft.fft(a, norm="ortho")


def _F2(a):
    return np.fft.ifft2(a, norm="ortho")


def _IF2(a):
    return np.fft.fft2(a, norm="ortho")


# ---- 1. Kuramoto-Sivashinsky, u_t = -u u_x - u_xx - u_xxxx ---------------

KSN = 2048
KSLL = 64.0 * np.pi        # domain length; large enough to be chaotic
KSDT = 0.01

ksxx = KSLL * np.arange(0, KSN) / KSN

# Wavenumbers in FFT order: 0, 1, ..., n/2-1, -n/2, ..., -1, scaled.
_j = np.arange(1, KSN + 1)
kskk = (2.0 * np.pi / KSLL) * np.where(_j <= KSN / 2, _j - 1, _j - 1 - KSN)

# Linear operator: -u_xx - u_xxxx becomes kappa^2 - kappa^4.  Stiff (kappa^4
# reaches 2.7e8), hence the implicit treatment below.
kslin = kskk ** 2 - kskk ** 4
ksden = 1.0 / (1.0 - KSDT * kslin)

ksu0 = np.cos(ksxx / 16.0) * (1.0 + np.sin(ksxx / 16.0))
ksh0 = _F(ksu0)


def ksstep(uh):
    """IMEX Euler: stiff linear part implicit, nonlinear part explicit.

    ``-u u_x`` is written as ``-1/2 d/dx (u^2)``, which costs one transform
    instead of two and is exactly conservative.
    """
    u = _IF(uh).real
    nl = (0.5j * kskk) * _F(u * u)
    return (uh + KSDT * nl) * ksden


def ksrun(m):
    return float((np.abs(nest(ksstep, ksh0, m)) ** 2).sum())


# ---- 2. 2D Navier-Stokes in vorticity form -------------------------------

NSN = 128
NSDT = 0.002
NSNU = 0.002               # viscosity

_j2 = np.arange(1, NSN + 1)
nsk1 = np.where(_j2 <= NSN / 2, _j2 - 1, _j2 - 1 - NSN).astype(float)
nskx = np.tile(nsk1, (NSN, 1))
nsky = np.tile(nsk1[:, None], (1, NSN))
nsk2 = nskx ** 2 + nsky ** 2

# The (0,0) mode is the mean, which the Poisson solve cannot determine; add 1
# there so the reciprocal is finite.  The mode is never used.
_e = np.zeros((NSN, NSN))
_e[0, 0] = 1.0
nsk2i = 1.0 / (nsk2 + _e)

# 2/3 dealiasing: a quadratic nonlinearity aliases the top third of the
# spectrum back onto the resolved modes, so those coefficients are zeroed.
nsdeal = (nsk1[:, None] ** 2 + nsk1[None, :] ** 2 <= (NSN / 3) ** 2).astype(float)

_ii = np.arange(1, NSN + 1)[:, None]
_jj = np.arange(1, NSN + 1)[None, :]
nsw0 = _F2(np.sin(2 * np.pi * _ii / NSN) * np.cos(4 * np.pi * _jj / NSN)
           + 0.4 * np.sin(6 * np.pi * (_ii + _jj) / NSN))


def nsstep(wh):
    """Stream function from vorticity by a Poisson solve, velocity from the
    stream function, then the advective term -- six transforms per step."""
    ph = wh * nsk2i                                # -Laplacian^-1 w
    u = _IF2((-1j * nsky) * ph).real               #  d psi / dy
    v = _IF2((1j * nskx) * ph).real                # -d psi / dx
    wx = _IF2((-1j * nskx) * wh).real
    wy = _IF2((-1j * nsky) * wh).real
    nl = _F2(u * wx + v * wy) * nsdeal
    return (wh - NSDT * nl) / (1.0 + NSDT * NSNU * nsk2)   # implicit viscosity


def nsrun(m):
    return float((np.abs(nest(nsstep, nsw0, m)) ** 2).sum())


# ---- 3. FFT Poisson solve ------------------------------------------------

PSN = 512
_pi = np.arange(1, PSN + 1)[:, None]
_pj = np.arange(1, PSN + 1)[None, :]
pssrc = np.sin(2 * np.pi * _pi / PSN) * np.cos(6 * np.pi * _pj / PSN)
_p1 = np.arange(1, PSN + 1)
psk1 = np.where(_p1 <= PSN / 2, _p1 - 1, _p1 - 1 - PSN).astype(float)
psk2 = psk1[:, None] ** 2 + psk1[None, :] ** 2
_pe = np.zeros((PSN, PSN))
_pe[0, 0] = 1.0
psk2i = 1.0 / (psk2 + _pe)


def psolve(m):
    """The source is varied per iteration so that no cache can answer twice."""
    r = 0.0
    k = 0
    while k < m:
        r = float((np.abs(_IF2(_F2(pssrc * (1.0 + 0.001 * k)) * psk2i)) ** 2).sum())
        k += 1
    return r


def main():
    print("Experiment 16 -- spectral PDE")
    print("")

    bench("Kuramoto-Sivashinsky, 2048 modes, 2000 steps", lambda: ksrun(2000))
    check("Kuramoto-Sivashinsky (50 steps)", ksrun(50))

    bench("2D Navier-Stokes vorticity, 128^2, 200 steps", lambda: nsrun(200))
    check("2D Navier-Stokes (20 steps)", nsrun(20))

    bench("FFT Poisson solve, 512^2, 30 solves", lambda: psolve(30))
    check("FFT Poisson solve (2 solves)", psolve(2))

    # Where a time step's cost actually is.  The transforms are not it.
    print("")
    print("-- one KS step, split --")
    kszh = ksh0
    bench("ifft (one transform)", lambda: _IF(kszh))
    kszu = _IF(kszh).real
    bench("fft (one transform)", lambda: _F(kszu * kszu))
    bench("one complex elementwise pass", lambda: kszh * ksden)
    bench("the whole step", lambda: ksstep(kszh))


if __name__ == "__main__":
    main()
