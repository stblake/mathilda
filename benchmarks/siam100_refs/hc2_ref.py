"""hc2 (holdout A, family C2) reference, method A: event-driven billiard in
mpmath (closed-form smaller root of |p + t v - c|^2 = r^2 over the lattice
points carrying a mirror whose distance to the ray is < r; v' = v - 2 (v.n) n).
Run at dps 150 and 200 (the agreement bounds the chaotic loss) and at dps 16
(the digits a perfect double computation keeps).  Method B: hc2_ref.wl
(independent WL implementation with NSolve hits at precision 150)."""
import mpmath as mp

F = mp.mpf
# k: (p0, v0 (unit, exact), r, T, mirror selector, quantity)
SEL = {"all": lambda i, j: True, "even": lambda i, j: (i + j) % 2 == 0,
       "oddcol": lambda i, j: i % 2 != 0}
PROBS = {
    1: ((F(1)/2, F(1)/10), (1, 0), F(1)/3, 1, "all", "norm"),
    2: ((F(1)/2, F(1)/10), (1, 0), F(1)/4, 4, "all", "norm"),
    3: ((F(1)/2, F(1)/2), (F(3)/5, F(4)/5), F(1)/3, 3, "all", "x"),
    4: ((F(0), F(1)/2), (F(5)/13, F(12)/13), F(2)/5, 3, "all", "norm"),
    5: ((F(1)/2, F(1)/10), (1, 0), F(1)/3, 2, "all", "norm"),
    6: ((F(1)/2, F(1)/2), (F(4)/5, F(3)/5), F(1)/2, 5, "even", "norm"),
    7: ((F(1)/2, F(1)/3), (F(8)/17, F(15)/17), F(1)/5, 8, "all", "y"),
    8: ((F(1)/2, F(1)/7), (-F(7)/25, F(24)/25), F(2)/5, 3, "all", "fromstart"),
    9: ((F(0), F(1)/2), (F(12)/13, F(5)/13), F(1)/3, 6, "oddcol", "norm"),
    10: ((F(1)/2, F(1)/10), (1, 0), F(1)/3, 5, "all", "norm"),
}


def run(k, dps):
    mp.mp.dps = dps
    (p0x, p0y), (v0x, v0y), r, T, sel, q = PROBS[k]
    # re-create the exact inputs at this precision
    p0x, p0y, r = mp.mpf(p0x), mp.mpf(p0y), mp.mpf(r)
    P0 = {1: (1, 2, 1, 10), 2: (1, 2, 1, 10), 3: (1, 2, 1, 2), 4: (0, 1, 1, 2), 5: (1, 2, 1, 10),
          6: (1, 2, 1, 2), 7: (1, 2, 1, 3), 8: (1, 2, 1, 7), 9: (0, 1, 1, 2), 10: (1, 2, 1, 10)}[k]
    V0 = {1: (1, 1, 0, 1), 2: (1, 1, 0, 1), 3: (3, 5, 4, 5), 4: (5, 13, 12, 13), 5: (1, 1, 0, 1),
          6: (4, 5, 3, 5), 7: (8, 17, 15, 17), 8: (-7, 25, 24, 25), 9: (12, 13, 5, 13), 10: (1, 1, 0, 1)}[k]
    R = {1: (1, 3), 2: (1, 4), 3: (1, 3), 4: (2, 5), 5: (1, 3), 6: (1, 2), 7: (1, 5), 8: (2, 5),
         9: (1, 3), 10: (1, 3)}[k]
    px, py = mp.mpf(P0[0]) / P0[1], mp.mpf(P0[2]) / P0[3]
    vx, vy = mp.mpf(V0[0]) / V0[1], mp.mpf(V0[2]) / V0[3]
    r = mp.mpf(R[0]) / R[1]
    sx, sy = px, py
    ok = SEL[sel]
    left = mp.mpf(T); last = None
    while True:
        best = None
        x0, x1 = sorted([px, px + left*vx]); y0, y1 = sorted([py, py + left*vy])
        for i in range(int(mp.floor(x0)) - 1, int(mp.ceil(x1)) + 2):
            for j in range(int(mp.floor(y0)) - 1, int(mp.ceil(y1)) + 2):
                if (i, j) == last or not ok(i, j):
                    continue
                dx, dy = i - px, j - py
                b = dx*vx + dy*vy
                if b <= 0:
                    continue
                disc = b*b - (dx*dx + dy*dy - r*r)
                if disc <= 0:
                    continue
                t = b - mp.sqrt(disc)
                if t > 0 and (best is None or t < best[0]):
                    best = (t, i, j)
        if best is None or best[0] >= left:
            px += left*vx; py += left*vy
            break
        t, i, j = best
        px += t*vx; py += t*vy; left -= t; last = (i, j)
        nx, ny = (px - i)/r, (py - j)/r
        d = vx*nx + vy*ny
        vx, vy = vx - 2*d*nx, vy - 2*d*ny
    return {"norm": mp.sqrt(px*px + py*py), "x": px, "y": py,
            "fromstart": mp.sqrt((px - sx)**2 + (py - sy)**2)}[q]


def nrefl(k):
    return None


if __name__ == "__main__":
    for k in sorted(PROBS):
        a = run(k, 150); b = run(k, 200); c = run(k, 16)
        mp.mp.dps = 200
        ag = -mp.log10(abs(a - b)/abs(b)) if a != b else 150
        d16 = -mp.log10(abs(c - b)/abs(b)) if c != b else 16
        d60 = run(k, 60); mp.mp.dps = 200
        d60 = -mp.log10(abs(d60 - b)/abs(b))
        print(k, mp.nstr(b, 70), "agree150/200: %.1f  dps16 keeps %.1f  dps60 keeps %.1f" % (float(ag), float(d16), float(d60)))
