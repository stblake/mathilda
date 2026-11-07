"""Shared reference helper for holdout-B families d11-d15 (experiments 246-250).

An arbitrary-precision Gragg-Bulirsch-Stoer integrator (modified midpoint +
polynomial extrapolation in h^2, adaptive macro-step).  It is used as the
SECOND, algorithmically independent method next to mpmath.odefun (a Taylor
series integrator) for the ODE references that have no closed form.

    bs_solve(f, t0, y0, t1, tol)   -> y(t1)   (f(t, y) -> list of mpf)

The macro-step is adapted so that the extrapolation tableau converges to tol
within KMAX columns; the error estimate is the difference of the last two
diagonal entries.  This is a reference tool only, never a benchmark case body.
"""
from mpmath import mp, mpf


def _midpoint(f, t, y, H, n):
    h = H / n
    z0 = list(y)
    fy = f(t, z0)
    z1 = [a + h * b for a, b in zip(z0, fy)]
    for m in range(1, n):
        fz = f(t + m * h, z1)
        z0, z1 = z1, [a + 2 * h * b for a, b in zip(z0, fz)]
    fz = f(t + H, z1)
    return [(a + b + h * c) / 2 for a, b, c in zip(z1, z0, fz)]


SEQ = [2 * (j + 1) for j in range(40)]          # 2, 4, 6, ... (Deuflhard)


def bs_step(f, t, y, H, tol, kmax):
    T = []
    for k in range(kmax):
        row = [_midpoint(f, t, y, H, SEQ[k])]
        for j in range(1, k + 1):
            r = (mpf(SEQ[k]) / SEQ[k - j]) ** 2 - 1
            row.append([a + (a - b) / r for a, b in zip(row[j - 1], T[k - 1][j - 1])])
        T.append(row)
        if k >= 3:
            err = max(abs(a - b) for a, b in zip(row[k], row[k - 1]))
            scale = max(1, max(abs(a) for a in row[k]))
            if err <= tol * scale:
                return row[k], k, err
    return None, kmax, err


def bs_solve(f, t0, y0, t1, tol=None, H=None, kmax=18):
    if tol is None:
        tol = mpf(10) ** (-(mp.dps - 8))
    t = mpf(t0)
    y = [mpf(v) for v in y0]
    t1 = mpf(t1)
    if H is None:
        H = (t1 - t) / 20
    while (t1 - t) * (1 if t1 > t0 else -1) > 0:
        if abs(H) > abs(t1 - t):
            H = t1 - t
        ynew, k, err = bs_step(f, t, y, H, tol, kmax)
        if ynew is None:
            H = H / 2
            continue
        t, y = t + H, ynew
        if k < kmax - 6:
            H = H * mpf(1.5)
        elif k > kmax - 3:
            H = H * mpf(0.7)
    return y
