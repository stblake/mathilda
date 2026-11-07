"""c2 reference (method 1): event-driven billiard in mpmath at high precision.
Photon from (1/2, 1/10) heading east at speed 1; mirrors of radius 1/3 at every
integer lattice point.  Each free flight is a straight line; the next hit is the
smaller root of |p + t v - c|^2 = r^2 in closed form over the lattice points
whose distance to the ray is < r; reflection v' = v - 2 (v.n) n.
Run at several precisions: the loss of digits (~38 at T = 10, book App. B)
is read off the agreement.  Prints |p(T)| for T in TS."""
import mpmath as mp, sys

def dist(T, dps):
    mp.mp.dps = dps
    r = mp.mpf(1)/3
    px, py = mp.mpf(1)/2, mp.mpf(1)/10
    vx, vy = mp.mpf(1), mp.mpf(0)
    left = mp.mpf(T); last = None
    while True:
        best = None
        # integer points near the segment p .. p + left v
        x0, x1 = sorted([px, px + left*vx]); y0, y1 = sorted([py, py + left*vy])
        for i in range(int(mp.floor(x0)) - 1, int(mp.ceil(x1)) + 2):
            for j in range(int(mp.floor(y0)) - 1, int(mp.ceil(y1)) + 2):
                if (i, j) == last: continue
                dx, dy = i - px, j - py
                b = dx*vx + dy*vy
                if b <= 0: continue
                disc = b*b - (dx*dx + dy*dy - r*r)
                if disc <= 0: continue
                t = b - mp.sqrt(disc)
                if t > 0 and (best is None or t < best[0]): best = (t, i, j)
        if best is None or best[0] >= left:
            px += left*vx; py += left*vy
            return mp.sqrt(px*px + py*py)
        t, i, j = best
        px += t*vx; py += t*vy; left -= t; last = (i, j)
        nx, ny = (px - i)/r, (py - j)/r
        d = vx*nx + vy*ny
        vx, vy = vx - 2*d*nx, vy - 2*d*ny

if __name__ == "__main__":
    for T in [1, 2, 3, 10]:
        a = dist(T, 150); b = dist(T, 200)
        print(T, mp.nstr(a, 70))
        print(T, mp.nstr(b, 70), "agree:", -int(mp.log10(abs(a-b)/b)) if a != b else ">150")
    for d in [16, 30, 40, 50, 60]:
        a = dist(10, d); mp.mp.dps = 200; t = dist(10, 200)
        print("dps", d, "T=10 digits:", mp.nstr(-mp.log10(abs(a-t)/t), 4))
    for T in [1, 1.5, 2, 3]:
        a = dist(T, 16); t = dist(T, 200)
        print("dps16 T=", T, "digits:", mp.nstr(-mp.log10(abs(a-t)/t), 4))
