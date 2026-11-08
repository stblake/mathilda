"""Experiment 227 -- SIAM 100-digit holdout B for Problem 2 (Python column).

The same statement as billiards.m: straight flight, the next hit as the
smallest positive real root of the intersection polynomial from the
library's polynomial root finder (numpy.roots in doubles, mpmath.polyroots
at 60 digits), reflection v' = v - 2 (v.n) n, the same candidate filter."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import numpy as np
import mpmath

require(["numpy:roots", "mpmath:polyroots"])

INF = None
T = {1: 6, 2: 5, 3: 4, 4: 12, 5: 6, 6: 20, 8: 100, 9: 3, 10: 25}
N7 = 16


def scene(k, F, sqrt):
    """(lattice, mirrors, p0, v0, T, reflections); mirror = (kind, a, b, c)"""
    h, s3 = F(1) / 2, sqrt(F(3))
    box = lambda x0, x1, y0, y1: [("line", F(1), F(0), F(x1)), ("line", F(-1), F(0), F(-x0)),
                                   ("line", F(0), F(1), F(y1)), ("line", F(0), F(-1), F(-y0))]
    sinai = box(-1, 1, -1, 1) + [("disk", F(0), F(0), h)]
    third, fifth = F(1) / 3, F(1) / 5
    return {
        1: (("sq", third), [], (h, F(1) / 10), (3 * fifth, 4 * fifth), T[1], INF),
        2: (("sq", F(1) / 4), [], (h, h), (4 * fifth, 3 * fifth), T[2], INF),
        3: (("tri", third), [], (h, F(1) / 5), (F(1), F(0)), T[3], INF),
        4: (None, sinai, (F(3) / 4, F(0)), (3 * fifth, 4 * fifth), T[4], INF),
        5: (None, box(0, 2, 0, 1) + [("disk", h, h, F(1) / 4), ("disk", F(3) / 2, h, F(1) / 4)],
            (F(1), F(1) / 10), (3 * fifth, 4 * fifth), T[5], INF),
        6: (None, [("disk", F(0), F(6) / 5, F(1)), ("disk", -3 * s3 / 5, F(-3) / 5, F(1)),
                   ("disk", 3 * s3 / 5, F(-3) / 5, F(1))],
            (F(0), F(0)), (F(12) / 13, F(-5) / 13), T[6], INF),
        7: (None, sinai, (-h, -F(3) / 4), (F(5) / 13, F(12) / 13), INF, N7),
        8: (None, box(0, 3, 0, 2), (F(1) / 3, F(1) / 7), (F(12) / 13, F(5) / 13), T[8], INF),
        9: (("sq", F(9) / 20), [], (h, h), (3 * fifth, 4 * fifth), T[9], INF),
        10: (None, sinai, (F(3) / 4, F(0)), (3 * fifth, 4 * fifth), T[10], INF),
    }[k]


def lattice(lat, p, q, F, sqrt, floor):
    kind, r = lat
    xs, ys = (p[0], q[0]), (p[1], q[1])
    if kind == "sq":
        return [(("disk", F(i), F(j), r), (i, j))
                for i in range(floor(min(xs)) - 1, floor(max(xs)) + 3)
                for j in range(floor(min(ys)) - 1, floor(max(ys)) + 3)]
    hh = sqrt(F(3)) / 2
    return [(("disk", i + F(j) / 2, j * hh, r), (i, j))
            for j in range(floor(min(ys) / hh) - 1, floor(max(ys) / hh) + 3)
            for i in range(floor(min(xs) - F(j) / 2) - 1, floor(max(xs) - F(j) / 2) + 3)]


def photon(k, num):
    F, roots, sqrt, floor = num
    lat, mir, (px, py), (vx, vy), Tk, nmax = scene(k, F, sqrt)
    left = F(Tk) if Tk is not None else None
    fin = [(m, ("w", i)) for i, m in enumerate(mir)]
    last, clock, nb = None, F(0), 0
    while True:
        span = F(50) if left is None else min(left, F(50))
        cands = fin + (lattice(lat, (px, py), (px + span * vx, py + span * vy), F, sqrt, floor) if lat else [])
        hits = []
        for (kind, a, b, c), key in cands:
            if key == last:
                continue
            if kind == "disk":
                dx, dy = a - px, b - py
                if dx * vx + dy * vy <= 0 or abs(vx * dy - vy * dx) >= c:
                    continue
                poly = [vx * vx + vy * vy, -2 * (dx * vx + dy * vy), dx * dx + dy * dy - c * c]
            else:
                if a * vx + b * vy == 0:
                    continue
                poly = [a * vx + b * vy, a * px + b * py - c]
            for t in roots(poly):
                if t > 0:
                    hits.append((t, kind, a, b, c, key))
        if not hits or (left is not None and min(hits)[0] >= left):
            px, py = px + left * vx, py + left * vy
            return sqrt(px * px + py * py)
        t, kind, a, b, c, key = min(hits)
        px, py = px + t * vx, py + t * vy
        clock += t; nb += 1; last = key
        if left is not None:
            left -= t
        nx, ny = ((px - a) / c, (py - b) / c) if kind == "disk" else (a, b)
        d = vx * nx + vy * ny
        vx, vy = vx - 2 * d * nx, vy - 2 * d * ny
        if nb == nmax:
            return clock


def np_real_roots(c):
    return [r.real for r in np.roots(c) if r.imag == 0]


def mp_real_roots(c):
    return [mpmath.re(r) for r in mpmath.polyroots(c, maxsteps=200, extraprec=60)
            if mpmath.im(r) == 0]


MACH = (float, np_real_roots, math.sqrt, lambda x: int(math.floor(x)))
MP = (mpmath.mpf, mp_real_roots, mpmath.sqrt, lambda x: int(mpmath.floor(x)))

for k in range(1, 11):
    vrun("bc2.%d machine" % k, lambda: photon(k, MACH), 10, 3, cap=60)
    vrun("bc2.%d wp50" % k, lambda: photon(k, MP), 48, cap=120, dps=60)
