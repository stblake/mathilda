"""Holdout B, family C2 (experiment 227) -- reference values, route A (mpmath).

Event-driven billiards: a unit-speed photon among fixed mirrors (disks, seen
from outside, and straight walls), specular reflection v' = v - 2 (v.n) n.
Next hit by the closed-form quadratic (disk) / linear (wall) root, exactly as
in siam100_refs/c2_ref.py.  Each problem is run at dps 150 and dps 250 (the
digits on which they agree are reported); and at dps 16 to measure how many
digits a double computation keeps.  Route B: siam100_refs/bc2_ref.gp (PARI).
"""
import sys
import mpmath as mp


def scene(k, F):
    """(lattice, finite obstacles, start p, direction v, T or None, nbounce or None)"""
    h = F(1) / 2
    s3 = mp.sqrt(3)
    box = lambda x0, x1, y0, y1: [("line", F(1), F(0), F(x1)), ("line", F(-1), F(0), -F(x0)),
                                   ("line", F(0), F(1), F(y1)), ("line", F(0), F(-1), -F(y0))]
    if k == 1: return (("sq", F(1) / 3), [], (h, F(1) / 10), (F(3) / 5, F(4) / 5), 6, None)
    if k == 2: return (("sq", F(1) / 4), [], (h, h), (F(4) / 5, F(3) / 5), 5, None)
    if k == 3: return (("tri", F(1) / 3), [], (h, F(1) / 5), (F(1), F(0)), 4, None)
    if k == 4: return (None, box(-1, 1, -1, 1) + [("disk", F(0), F(0), h)],
                       (F(3) / 4, F(0)), (F(3) / 5, F(4) / 5), 12, None)
    if k == 5: return (None, box(0, 2, 0, 1) + [("disk", h, h, F(1) / 4), ("disk", F(3) / 2, h, F(1) / 4)],
                       (F(1), F(1) / 10), (F(3) / 5, F(4) / 5), 6, None)
    if k == 6: return (None, [("disk", F(0), F(6) / 5, F(1)), ("disk", -3 * s3 / 5, -F(3) / 5, F(1)),
                              ("disk", 3 * s3 / 5, -F(3) / 5, F(1))],
                       (F(0), F(0)), (F(12) / 13, -F(5) / 13), 20, None)
    if k == 7: return (None, box(-1, 1, -1, 1) + [("disk", F(0), F(0), h)],
                       (-h, -F(3) / 4), (F(5) / 13, F(12) / 13), None, 16)
    if k == 8: return (None, box(0, 3, 0, 2), (F(1) / 3, F(1) / 7), (F(12) / 13, F(5) / 13), 100, None)
    if k == 9: return (("sq", F(9) / 20), [], (h, h), (F(3) / 5, F(4) / 5), 3, None)
    if k == 10: return (None, box(-1, 1, -1, 1) + [("disk", F(0), F(0), h)],
                        (F(3) / 4, F(0)), (F(3) / 5, F(4) / 5), 25, None)


def lattice_disks(lat, p, q):
    kind, r = lat
    xs, ys = (p[0], q[0]), (p[1], q[1])
    out = []
    if kind == "sq":
        for i in range(int(mp.floor(min(xs))) - 1, int(mp.floor(max(xs))) + 3):
            for j in range(int(mp.floor(min(ys))) - 1, int(mp.floor(max(ys))) + 3):
                out.append(("disk", mp.mpf(i), mp.mpf(j), r, (i, j)))
    else:                                    # triangular: (i + j/2, j sqrt3/2)
        hh = mp.sqrt(3) / 2
        for j in range(int(mp.floor(min(ys) / hh)) - 1, int(mp.floor(max(ys) / hh)) + 3):
            for i in range(int(mp.floor(min(xs) - j / mp.mpf(2))) - 1, int(mp.floor(max(xs) - j / mp.mpf(2))) + 3):
                out.append(("disk", i + j / mp.mpf(2), j * hh, r, (i, j)))
    return out


def fly(k):
    F = mp.mpf
    lat, obs, p, v, T, nb = scene(k, F)
    px, py = p; vx, vy = v
    left = F(T) if T is not None else F(10) ** 9
    last, clock, bounces = None, F(0), 0
    obs = [o + (("w", idx),) for idx, o in enumerate(obs)]
    while True:
        q = (px + min(left, 50) * vx, py + min(left, 50) * vy)
        cands = list(obs) + (lattice_disks(lat, (px, py), q) if lat else [])
        best = None
        for o in cands:
            if o[-1] == last:
                continue
            if o[0] == "disk":
                _, cx, cy, r, key = o
                dx, dy = cx - px, cy - py
                if dx * vx + dy * vy <= 0 or abs(vx * dy - vy * dx) >= r:
                    continue
                b = -2 * (dx * vx + dy * vy); c = dx * dx + dy * dy - r * r
                a = vx * vx + vy * vy                   # not assumed 1: keeps hits on the circle
                disc = b * b - 4 * a * c
                if disc < 0:
                    continue
                t = (-b - mp.sqrt(disc)) / (2 * a)
            else:
                _, nx, ny, cc, key = o
                den = nx * vx + ny * vy
                if den == 0:
                    continue
                t = (cc - nx * px - ny * py) / den
            if t > 0 and (best is None or t < best[0]):
                best = (t, o)
        if best is None or (T is not None and best[0] >= left):
            px, py = px + left * vx, py + left * vy
            return mp.sqrt(px * px + py * py)
        t, o = best
        px, py = px + t * vx, py + t * vy
        left -= t; clock += t; bounces += 1; last = o[-1]
        if o[0] == "disk":
            nx, ny = (px - o[1]) / o[3], (py - o[2]) / o[3]
        else:
            nx, ny = o[1], o[2]
        d = vx * nx + vy * ny
        vx, vy = vx - 2 * d * nx, vy - 2 * d * ny
        if nb is not None and bounces == nb:
            return clock


if __name__ == "__main__":
    ks = [int(a) for a in sys.argv[1:]] or range(1, 11)
    for k in ks:
        mp.mp.dps = 250; hi = fly(k)
        mp.mp.dps = 150; md = fly(k)
        mp.mp.dps = 16; lo = fly(k)
        mp.mp.dps = 60; w = fly(k)
        mp.mp.dps = 260
        ag = lambda a: float(-mp.log10(abs(a - hi) / abs(hi))) if a != hi else 99
        print("bc2.%d %s\n   agree(150 vs 250) %.1f  dps16 keeps %.1f  dps60 keeps %.1f" %
              (k, mp.nstr(hi, 80), ag(md), ag(lo), ag(w)), flush=True)
