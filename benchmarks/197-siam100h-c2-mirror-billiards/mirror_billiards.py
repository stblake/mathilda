"""Experiment 197 -- SIAM 100-digit holdout A, family hc2 (Python column).

The statement of mirror_billiards.m (C2's, experiment 167, with the start,
heading, mirror radius, mirror set, time and reported quantity varied): the
next hit as the smallest positive real root of |p + t v - c|^2 = r^2 from the
library's polynomial root finder (numpy.roots in doubles, mpmath.polyroots at
60 digits), reflection v' = v - 2 (v.n) n, over the same candidate filter."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
from fractions import Fraction as Q
import math
import numpy as np
import mpmath

require(["numpy:roots", "mpmath:polyroots"])

SEL = {"all": lambda i, j: True, "even": lambda i, j: (i + j) % 2 == 0,
       "oddcol": lambda i, j: i % 2 != 0}
HC2 = [((Q(1, 2), Q(1, 10)), (Q(1), Q(0)), Q(1, 3), 1, "all", "norm"),
       ((Q(1, 2), Q(1, 10)), (Q(1), Q(0)), Q(1, 4), 4, "all", "norm"),
       ((Q(1, 2), Q(1, 2)), (Q(3, 5), Q(4, 5)), Q(1, 3), 3, "all", "x"),
       ((Q(0), Q(1, 2)), (Q(5, 13), Q(12, 13)), Q(2, 5), 3, "all", "norm"),
       ((Q(1, 2), Q(1, 10)), (Q(1), Q(0)), Q(1, 3), 2, "all", "norm"),
       ((Q(1, 2), Q(1, 2)), (Q(4, 5), Q(3, 5)), Q(1, 2), 5, "even", "norm"),
       ((Q(1, 2), Q(1, 3)), (Q(8, 17), Q(15, 17)), Q(1, 5), 8, "all", "y"),
       ((Q(1, 2), Q(1, 7)), (Q(-7, 25), Q(24, 25)), Q(2, 5), 3, "all", "fromstart"),
       ((Q(0), Q(1, 2)), (Q(12, 13), Q(5, 13)), Q(1, 3), 6, "oddcol", "norm"),
       ((Q(1, 2), Q(1, 10)), (Q(1), Q(0)), Q(1, 3), 5, "all", "norm")]
NEED_M = [10, 10, 10, 10, 10, 10, 6, 10, 10, 10]


def photon(prob, num):
    """num = (mk(Fraction), real roots of a quadratic [a, b, c], sqrt, floor)"""
    (p0x, p0y), (v0x, v0y), r0, T, s, qty = prob
    mk, roots, sqrt, floor = num
    ok = SEL[s]
    r = mk(r0)
    px, py, vx, vy, left, last = mk(p0x), mk(p0y), mk(v0x), mk(v0y), mk(T), None
    sx, sy = px, py
    while True:
        xs, ys = (px, px + left * vx), (py, py + left * vy)
        hits = []
        for i in range(floor(min(xs)) - 1, floor(max(xs)) + 3):
            for j in range(floor(min(ys)) - 1, floor(max(ys)) + 3):
                dx, dy = i - px, j - py
                if (not ok(i, j) or (i, j) == last or dx * vx + dy * vy <= 0
                        or abs(vx * dy - vy * dx) >= r):
                    continue
                for t in roots([vx * vx + vy * vy, -2 * (dx * vx + dy * vy),
                                dx * dx + dy * dy - r * r]):
                    if t > 0:
                        hits.append((t, i, j))
        if not hits or min(hits)[0] >= left:
            px, py = px + left * vx, py + left * vy
            break
        t, i, j = min(hits)
        px, py, left, last = px + t * vx, py + t * vy, left - t, (i, j)
        nx, ny = (px - i) / r, (py - j) / r
        d = vx * nx + vy * ny
        vx, vy = vx - 2 * d * nx, vy - 2 * d * ny
    return {"norm": lambda: sqrt(px * px + py * py), "x": lambda: px, "y": lambda: py,
            "fromstart": lambda: sqrt((px - sx) ** 2 + (py - sy) ** 2)}[qty]()


def np_real_roots(c):
    return [r.real for r in np.roots(c) if r.imag == 0]


def mp_real_roots(c):
    return [mpmath.re(r) for r in mpmath.polyroots(c, maxsteps=200, extraprec=60)
            if mpmath.im(r) == 0]


MACH = (lambda q: q.numerator / q.denominator if isinstance(q, Q) else float(q),
        np_real_roots, math.sqrt, lambda x: int(math.floor(x)))
MP = (lambda q: mpmath.mpf(q.numerator) / q.denominator if isinstance(q, Q) else mpmath.mpf(q),
      mp_real_roots, mpmath.sqrt, lambda x: int(mpmath.floor(x)))

for k in range(10):
    vrun("hc2.%d machine" % (k + 1), lambda: photon(HC2[k], MACH), NEED_M[k], cap=60)
    vrun("hc2.%d wp50" % (k + 1), lambda: photon(HC2[k], MP), 48, cap=60, dps=60)
