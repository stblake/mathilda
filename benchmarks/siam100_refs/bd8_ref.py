# Holdout B, family D8 (rational integrands over R^N).
# Method A (this file, mpmath dps 70): spherical coordinates; the radial
# integral Int_0^inf r^(N-1)/P(r w) dr in closed form by residues
# (t = r^2, roots of the polynomial A(t) = sum_j h_j(w^2) t^j by polyroots),
# the angular integral over the positive orthant by tanh-sinh (N = 2) or
# Gauss-Legendre in cos(theta) x trapezoid in phi (N = 3; run at two
# resolutions, --spec M NT, to confirm convergence).
# Method B: bd8_ref.gp (PARI/gp: Cartesian iterated integrals, the innermost
# x-integral by residues, the outer ones by intnum) -- and closed forms for
# bd8.1, bd8.7-bd8.10.
import mpmath as mp
mp.mp.dps = 70
from itertools import combinations_with_replacement as cwr

def hj(u, j):
    """complete homogeneous symmetric polynomial h_j(u)"""
    if j == 0:
        return mp.mpf(1)
    s = mp.mpf(0)
    for c in cwr(range(len(u)), j):
        p = mp.mpf(1)
        for i in c:
            p *= u[i]
        s += p
    return s

def radial(coeffs, N):
    """Int_0^inf r^(N-1) / sum_j coeffs[j] r^(2j) dr  (= 1/2 Int_0^inf t^(N/2-1)/A(t) dt)"""
    A = coeffs[::-1]                       # highest first for polyroots
    dA = [c * (len(A) - 1 - i) for i, c in enumerate(A[:-1])]
    roots = mp.polyroots(A, maxsteps=200, extraprec=200)
    s = mp.mpf(N) / 2
    tot = 0
    for t in roots:
        d = mp.polyval(dA, t)
        if N % 2:
            tot += (-t) ** (s - 1) / d
        else:
            tot += t ** (s - 1) * mp.log(-t) / d
    if N % 2:
        return mp.re(mp.pi / mp.sin(mp.pi * s) * tot) / 2
    return mp.re(-tot) / 2

def sob2(m, N):
    """c_m(R^N)^2 * (2 pi)^N  = Int_{R^N} 1/sum_{|k|<=m} x^(2k)"""
    if N == 1:
        return 2 * radial([1] * (m + 1), 1)
    if N == 2:
        f = lambda th: radial([hj([mp.cos(th) ** 2, mp.sin(th) ** 2], j) for j in range(m + 1)], 2)
        return 4 * mp.quad(f, [0, mp.pi / 4, mp.pi / 2])
    if N == 3:
        # octant, t = cos(theta) in [0, 1] by Gauss-Legendre (the integrand is
        # a function of t^2: analytic), phi by the trapezoid rule over its
        # period pi (u1, u2 depend on cos 2 phi): both exponentially convergent.
        M, NT = SPEC3
        ts, ws = gl(NT)
        tot = 0
        for t, wt in zip(ts, ws):
            for i in range(M):
                ph = mp.pi * i / M
                s2 = 1 - t * t
                u = [s2 * mp.cos(ph) ** 2, s2 * mp.sin(ph) ** 2, t * t]
                tot += wt * radial([hj(u, j) for j in range(m + 1)], 3)
        return 8 * tot * (mp.pi / M) / 2

SPEC3 = (120, 80)
def gl(n):
    """Gauss-Legendre nodes/weights on [0, 1]"""
    xs, ws = [], []
    for k in range(1, n + 1):
        x = mp.cos(mp.pi * (k - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp; x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        xs.append((1 + x) / 2); ws.append(1 / ((1 - x * x) * dp * dp))
    return xs, ws

def c(m, N):
    return mp.sqrt(sob2(m, N) / (2 * mp.pi) ** N)

import sys
if len(sys.argv) > 2 and sys.argv[1] == '--spec':
    SPEC3 = (int(sys.argv[2]), int(sys.argv[3])); sys.argv = sys.argv[:1] + sys.argv[4:]
which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
for k in which:
    if k == "1": print("bd8.1", mp.nstr(c(2, 1), 66), " closed", mp.nstr(1 / mp.sqrt(2 * mp.sqrt(3)) , 66))
    if k == "2": print("bd8.2", mp.nstr(c(6, 1), 66))
    if k == "3": print("bd8.3", mp.nstr(c(3, 2), 66))
    if k == "4": print("bd8.4", mp.nstr(c(4, 2), 66))
    if k == "5": print("bd8.5", mp.nstr(c(3, 3), 66))
    if k == "6": print("bd8.6", mp.nstr(c(2, 3), 66))
    if k == "7":   # Int_{R^2} 1/(1 + x^4 + y^4) = pi K(1/2)
        f = lambda th: radial([1, 0, mp.cos(th) ** 4 + mp.sin(th) ** 4], 2)
        print("bd8.7", mp.nstr(4 * mp.quad(f, [0, mp.pi / 4, mp.pi / 2]), 66), " closed", mp.nstr(mp.pi * mp.ellipk(mp.mpf(1) / 2), 66))
    if k == "8":   # Int_{R^3} (1 + |x|^2)^-3 = pi^2/4
        print("bd8.8", mp.nstr(4 * mp.pi * mp.quad(lambda r: r ** 2 / (1 + r ** 2) ** 3, [0, mp.inf]), 66), " closed", mp.nstr(mp.pi ** 2 / 4, 66))
    if k == "9":   # Int_{R^4} (1 + |x|^2)^-3 = pi^2/2
        print("bd8.9", mp.nstr(2 * mp.pi ** 2 * mp.quad(lambda r: r ** 3 / (1 + r ** 2) ** 3, [0, mp.inf]), 66), " closed", mp.nstr(mp.pi ** 2 / 2, 66))
    if k == "10":  # Int_{R^2} 1/(1 + x^4 + x^2 y^2 + y^4) = pi K(1/4)
        f = lambda th: radial([1, 0, mp.cos(th) ** 4 + (mp.cos(th) * mp.sin(th)) ** 2 + mp.sin(th) ** 4], 2)
        print("bd8.10", mp.nstr(4 * mp.quad(f, [0, mp.pi / 4, mp.pi / 2]), 66), " closed", mp.nstr(mp.pi * mp.ellipk(mp.mpf(1) / 4), 66))
    sys.stdout.flush()
