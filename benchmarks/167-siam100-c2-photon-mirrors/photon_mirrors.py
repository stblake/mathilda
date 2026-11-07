"""Experiment 167 -- SIAM 100-Digit Challenge, Problem 2 (Python column).

The same statement as photon_mirrors.m: straight-line flight, the next hit as
the smallest positive real root of |p + t v - c|^2 = 1/9 from the library's
polynomial root finder (numpy.roots in doubles, mpmath.polyroots at 60
digits), reflection v' = v - 2 (v.n) n, over the same candidate filter."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import numpy as np
import mpmath

require(["numpy:roots", "mpmath:polyroots"])


def photon(T, num):
    """num = (scalar(x), real roots of a quadratic [a, b, c], sqrt, floor)"""
    mk, roots, sqrt, floor = num
    third = mk(1) / 3
    px, py, vx, vy, left, last = mk(1) / 2, mk(1) / 10, mk(1), mk(0), mk(T), None
    while True:
        xs, ys = (px, px + left * vx), (py, py + left * vy)
        hits = []
        for i in range(floor(min(xs)) - 1, floor(max(xs)) + 3):
            for j in range(floor(min(ys)) - 1, floor(max(ys)) + 3):
                dx, dy = i - px, j - py
                if (i, j) == last or dx * vx + dy * vy <= 0 or abs(vx * dy - vy * dx) >= third:
                    continue
                # |p + t v - c|^2 - 1/9 = (v.v) t^2 - 2 (v.d) t + (d.d - 1/9)
                for t in roots([vx * vx + vy * vy, -2 * (dx * vx + dy * vy),
                                dx * dx + dy * dy - third * third]):
                    if t > 0:
                        hits.append((t, i, j))
        if not hits or min(hits)[0] >= left:
            px, py = px + left * vx, py + left * vy
            return sqrt(px * px + py * py)
        t, i, j = min(hits)
        px, py, left, last = px + t * vx, py + t * vy, left - t, (i, j)
        nx, ny = 3 * (px - i), 3 * (py - j)
        d = vx * nx + vy * ny
        vx, vy = vx - 2 * d * nx, vy - 2 * d * ny


def np_real_roots(c):
    return [r.real for r in np.roots(c) if r.imag == 0]


def mp_real_roots(c):
    return [mpmath.re(r) for r in mpmath.polyroots(c, maxsteps=200, extraprec=60)
            if mpmath.im(r) == 0]


MACH = (float, np_real_roots, math.sqrt, lambda x: int(math.floor(x)))
MP = (mpmath.mpf, mp_real_roots, mpmath.sqrt, lambda x: int(mpmath.floor(x)))

vrun("c2.machine", lambda: photon(10, MACH), 3, 3)
vrun("c2.machine-t3", lambda: photon(3, MACH), 10, 3)
vrun("c2.wp50", lambda: photon(10, MP), 45, dps=60)
