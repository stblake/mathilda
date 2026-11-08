"""References A and B for holdout-B family D16 (exp 251): the heat equation
from rough initial data; the time t* at which a functional of u(., t) takes a
given value.  u_t = u_xx on [0,1] (bd16.6: [0,pi]; bd16.9: the unit square).

  bd16.1  Dirichlet, u0 = 1:                          u(1/2, t) = 1/2
  bd16.2  Dirichlet, u0 = 1 - |2x - 1| (hat):          Int_0^1 u dx = 1/4
  bd16.3  Dirichlet, u0 = Haar psi (1 on [0,1/2), -1 on [1/2,1)): arc length of the graph = 3
  bd16.4  Dirichlet, u0 = x:                          arc length = 2
  bd16.5  Neumann both ends, u0 = [x > 1/3]:           u(0, t) = 1/3
  bd16.6  [0,pi], Dirichlet, u0 = x (pi - x):         Int_0^pi u^2 dx = 1
  bd16.7  u_t = u_xx + 1, Dirichlet, u0 = 0:          u(1/2, t) = 1/10
  bd16.8  u(0) = 0, u_x(1) + u(1) = 0 (Robin), u0 = 1: u(1, t) = 1/4
  bd16.9  unit square, Dirichlet, u0 = 1:              u(1/2, 1/2, t) = 1/2
  bd16.10 Dirichlet, u0 = Haar phi(3x - 1) = [1/3 <= x < 2/3]: arc length = 28/10

Method A: eigenfunction (Fourier) series with closed-form coefficients
  (bd16.8: Robin eigenvalues tan(mu) = -mu by findroot), the functional by
  term-wise closed forms (point values, integrals) or tanh-sinh quadrature of
  sqrt(1 + u_x^2) (arc length), the time by findroot.  bd16.9 = (bd16.1's
  u(1/2,t))^2 by separation of variables (u(x,y,t) = U(x,t) U(y,t)).
Method B: the method of images -- u = heat kernel * (odd / even periodic
  extension of u0) as a sum of erf / Gaussian terms for piecewise-linear data
  (bd16.7: tanh-sinh quadrature of kernel * data); bd16.6 (t* ~ 1.2, where the
  kernel is wider than the interval) instead by sine coefficients computed by
  quadrature and Int u^2 by quadrature of u^2 (no closed form, no Parseval);
  bd16.8 instead by
  numerical inversion of the closed-form Laplace transform (mpmath
  invertlaplace, Talbot) -- no eigenvalues.
Both at 70 digits; each prints t* and the number of agreeing digits.
Usage: python3.11 bd16_ref.py [labels...]
"""
import sys
import mpmath as mp

mp.mp.dps = 70
pi = mp.pi
half = mp.mpf(1) / 2
TOL = mp.mpf(10) ** -66


# ----------------------------------------------------------------- method A
def sine_series(b, lam, x, t, d=0, tol=TOL, nmin=5):
    """sum_n b(n) e^{-lam(n) t} sin^(d)(sqrt(lam(n)) x)  (truncated by size)"""
    s = mp.mpf(0); n = 1; small = 0
    while True:
        k = mp.sqrt(lam(n)); bn = b(n)
        term = bn * mp.exp(-lam(n) * t) * (mp.sin(k * x) if d == 0 else k * mp.cos(k * x))
        s += term
        if n > nmin and abs(bn * mp.exp(-lam(n) * t)) * (k if d else 1) < tol:
            small += 1
            if small > 3:
                return s
        else:
            small = 0
        n += 1


def arclen(ux, pts):
    return mp.quad(lambda x: mp.sqrt(1 + ux(x) ** 2), pts)


def cos_series_d1(b, t, tol=TOL):
    """u_x(x) = sum_n b(n) (n pi) e^{-(n pi)^2 t} cos(n pi x) on [0,1]: the coefficients for
    this t once, then each point by Clenshaw's recurrence (no per-term cosine)."""
    a = []; n = 1
    while True:
        an = b(n) * n * pi * mp.exp(-(n * pi) ** 2 * t)
        a.append(an)
        if n > 5 and mp.exp(-(n * pi) ** 2 * t) * n < tol:
            break
        n += 1

    def f(x):
        c = mp.cos(pi * x); b1 = b2 = mp.mpf(0)
        for an in reversed(a):          # sum_{n>=1} a_n cos(n th) = b1 cos th - b2, Clenshaw
            b1, b2 = an + 2 * c * b1 - b2, b1
        return b1 * c - b2
    return f


# Dirichlet sine coefficients on [0,1] (lam = (n pi)^2)
lamD = lambda n: (n * pi) ** 2
bBox = lambda n: 0 if n % 2 == 0 else 4 / (n * pi)
bHat = lambda n: 8 * mp.sin(n * pi / 2) / (n * pi) ** 2
bHaar = lambda n: 2 * (1 - 2 * mp.cos(n * pi / 2) + mp.cos(n * pi)) / (n * pi)
bX = lambda n: 2 * (-1) ** (n + 1) / (n * pi)
bBox3 = lambda n: 2 * (mp.cos(n * pi / 3) - mp.cos(2 * n * pi / 3)) / (n * pi)

A = {}
A["bd16.1"] = (lambda t: sine_series(bBox, lamD, half, t) - half, 0.1)
def _int_series(b, t):          # Int_0^1 sum b_n e^{-(n pi)^2 t} sin(n pi x) dx
    s = 0; n = 1
    while True:
        term = b(n) * mp.exp(-(n * pi) ** 2 * t) * (1 - mp.cos(n * pi)) / (n * pi)
        s += term
        if n > 5 and mp.exp(-(n * pi) ** 2 * t) / n < TOL:     # (b_n vanishes for even n: bound by size)
            return s
        n += 1


A["bd16.2"] = (lambda t: _int_series(bHat, t) - mp.mpf(1) / 4, 0.03)
A["bd16.3"] = (lambda t: arclen(cos_series_d1(bHaar, t), [0, half, 1]) - 3, 0.0155)
A["bd16.4"] = (lambda t: arclen(cos_series_d1(bX, t), [0, half, 1]) - 2, 0.005)


def _neu5(t):      # cosine series of [x > 1/3]: a0 = 2/3, a_n = -2 sin(n pi/3)/(n pi)
    s = mp.mpf(2) / 3; n = 1
    while True:
        an = -2 * mp.sin(n * pi / 3) / (n * pi)
        s += an * mp.exp(-(n * pi) ** 2 * t)
        if n > 5 and mp.exp(-(n * pi) ** 2 * t) / n < TOL:
            return s
        n += 1


A["bd16.5"] = (lambda t: _neu5(t) - mp.mpf(1) / 3, 0.02)


def _en6(t):       # x(pi-x) on [0,pi]: b_n = 8/(pi n^3) (odd n); Int u^2 = (pi/2) sum b_n^2 e^{-2 n^2 t}
    s = 0; n = 1
    while True:
        s += (8 / (pi * n ** 3)) ** 2 * mp.exp(-2 * n * n * t)
        if n > 5 and mp.exp(-2 * n * n * t) / n ** 6 < TOL:
            return pi / 2 * s
        n += 2


A["bd16.6"] = (lambda t: _en6(t) - 1, 1.0)
# bd16.7: u = x(1-x)/2 - v, v heat with v0 = x(1-x)/2, b_n = 4/(n pi)^3 (odd n)
A["bd16.7"] = (lambda t: mp.mpf(1) / 8 - sine_series(lambda n: 0 if n % 2 == 0 else 4 / (n * pi) ** 3, lamD, half, t)
               - mp.mpf(1) / 10, 0.1)


def robin_mus(N):     # roots of sin(mu) + mu cos(mu) = 0 in ((k - 1/2) pi, k pi)
    return [mp.findroot(lambda m: mp.sin(m) + m * mp.cos(m), ((k - half) * pi + mp.mpf(10) ** -20, k * pi),
                        solver="anderson") for k in range(1, N + 1)]


_MU = []


def _robin8(t):
    global _MU
    s = 0; k = 0
    while True:
        if k >= len(_MU):
            _MU += robin_mus(len(_MU) + 40)[len(_MU):]
        m = _MU[k]
        bn = (1 - mp.cos(m)) / m / (half - mp.sin(2 * m) / (4 * m))
        term = bn * mp.exp(-m * m * t) * mp.sin(m)
        s += term
        if k > 5 and abs(bn) * mp.exp(-m * m * t) < TOL:
            return s
        k += 1


A["bd16.8"] = (lambda t: _robin8(t) - mp.mpf(1) / 4, 0.3)
A["bd16.9"] = (lambda t: sine_series(bBox, lamD, half, t) ** 2 - half, 0.05)
A["bd16.10"] = (lambda t: arclen(cos_series_d1(bBox3, t), [0, mp.mpf(1) / 3, half, mp.mpf(2) / 3, 1])
                - mp.mpf(28) / 10, 2e-4)


# ----------------------------------------------------------------- method B
def G(s, t):
    return mp.exp(-s * s / (4 * t)) / mp.sqrt(4 * pi * t)


def extend(pieces, L, odd, K):
    """pieces (a, b, alpha, beta): u0 = alpha + beta y on (a,b) in [0,L]; return the
    pieces of the odd/even 2L-periodic extension on [-2LK, 2LK]."""
    out = []
    for k in range(-K, K + 1):
        sh = 2 * L * k
        for a, b, al, be in pieces:
            out.append((a + sh, b + sh, al - be * sh, be))                 # y -> y - sh
            # reflection y -> -y: value s * (al + be(-y)) on (-b, -a)
            sg = -1 if odd else 1
            out.append((-b + sh, -a + sh, sg * (al + be * sh), -sg * be))   # u0(-(y - sh))
    return out


def img_u(ext, x, t):
    s = 0; r = 2 * mp.sqrt(t)
    for a, b, al, be in ext:
        E = (mp.erf((x - a) / r) - mp.erf((x - b) / r)) / 2
        s += al * E + be * (x * E + mp.sqrt(t / pi) * (mp.exp(-(a - x) ** 2 / (4 * t)) - mp.exp(-(b - x) ** 2 / (4 * t))))
    return s


def near(ext, t, lo=0, hi=1):
    """the image pieces that matter on [lo, hi] at time t: farther than 30 sqrt(t) (where the
    kernel is < e^-225) they contribute nothing at 70 digits"""
    w = 30 * mp.sqrt(t)
    return [p for p in ext if p[1] > lo - w and p[0] < hi + w]


def img_ux(ext, x, t):
    s = 0; r = 2 * mp.sqrt(t)
    for a, b, al, be in near(ext, t):
        s += (al + be * a) * G(x - a, t) - (al + be * b) * G(x - b, t)
        if be:
            s += be * (mp.erf((x - a) / r) - mp.erf((x - b) / r)) / 2
    return s


def img_int(ext, t, lo, hi):    # Int_lo^hi u dx by quadrature of the erf sum
    return mp.quad(lambda x: img_u(ext, x, t), [lo, (lo + hi) / 2, hi])


one = mp.mpf(1)
P = {"box": [(0, one, one, 0)], "hat": [(0, half, 0, 2 * one), (half, one, 2 * one, -2 * one)],
     "haar": [(0, half, one, 0), (half, one, -one, 0)], "x": [(0, one, 0, one)],
     "step": [(one / 3, one, one, 0)], "box3": [(one / 3, 2 * one / 3, one, 0)]}
EXT = {k: extend(v, one, k != "step", 6) for k, v in P.items()}


def quad_images(f0, L, odd, x, t, K=6):
    """u(x,t) = sum_k Int_0^L [G(x-y-2Lk) -+ G(x+y-2Lk)] f0(y) dy  (tanh-sinh)"""
    sg = -1 if odd else 1
    tot = 0
    for k in range(-K, K + 1):
        c1, c2 = x - 2 * L * k, -(x - 2 * L * k)
        pts = sorted(set([mp.mpf(0), L] + [p for p in (c1, c2) if 0 < p < L]))
        tot += mp.quad(lambda y: (G(x - y - 2 * L * k, t) + sg * G(x + y - 2 * L * k, t)) * f0(y), pts)
    return tot


B = {}
B["bd16.1"] = lambda t: img_u(EXT["box"], half, t) - half
B["bd16.2"] = lambda t: img_int(EXT["hat"], t, 0, one) - one / 4
B["bd16.3"] = lambda t: arclen(lambda x: img_ux(EXT["haar"], x, t), [0, half, 1]) - 3
B["bd16.4"] = lambda t: arclen(lambda x: img_ux(EXT["x"], x, t), [0, half, 1]) - 2
B["bd16.5"] = lambda t: img_u(EXT["step"], mp.mpf(0), t) - one / 3
# bd16.6 at t* ~ 1.2 the kernel is wider than the interval (images converge slowly), so B
# is instead: sine coefficients by QUADRATURE (not the closed form), u summed on the fly, and
# Int u^2 by quadrature of u^2 (not Parseval).
_B6 = {}


def _b6(n):
    if n not in _B6:
        _B6[n] = 2 / pi * mp.quad(lambda y: y * (pi - y) * mp.sin(n * y), [0, pi / 2, pi])
    return _B6[n]


def _u6(x, t):
    s = 0; n = 1
    while True:
        s += _b6(n) * mp.exp(-n * n * t) * mp.sin(n * x)
        if n > 5 and mp.exp(-n * n * t) < TOL:
            return s
        n += 1


B["bd16.6"] = lambda t: mp.quad(lambda x: _u6(x, t) ** 2, [0, pi / 2, pi]) - 1
B["bd16.7"] = lambda t: one / 8 - quad_images(lambda y: y * (1 - y) / 2, one, True, half, t) - one / 10


def lap8(t):
    # u_hat(1, s) for u(0) = 0, u_x(1) + u(1) = 0, u0 = 1:  u_hat = 1/s + A sinh(q x) - cosh(q x)/s, q = sqrt s
    def F(s):
        q = mp.sqrt(s)
        Acoef = (q * mp.sinh(q) / s + mp.cosh(q) / s - 1 / s) / (q * mp.cosh(q) + mp.sinh(q))
        return 1 / s + Acoef * mp.sinh(q) - mp.cosh(q) / s
    return mp.invertlaplace(F, t, method="talbot")


B["bd16.8"] = lambda t: lap8(t) - one / 4
B["bd16.9"] = lambda t: img_u(EXT["box"], half, t) ** 2 - half
B["bd16.10"] = lambda t: arclen(lambda x: img_ux(EXT["box3"], x, t), [0, one / 3, half, 2 * one / 3, 1]) - mp.mpf(28) / 10


def solve(f, t0):
    # log-time secant: robust for the small-t arc-length problems
    s = mp.findroot(lambda s: f(mp.exp(s)), (mp.log(t0), mp.log(t0 * 1.1)), solver="secant", tol=mp.mpf(10) ** -64, maxsteps=200)
    return mp.exp(s)


if __name__ == "__main__":
    labels = sys.argv[1:] or ["bd16.%d" % k for k in range(1, 11)]
    for lab in labels:
        fA, t0 = A[lab]
        ta = solve(fA, t0)
        print(lab, "A", mp.nstr(ta, 66), flush=True)
        tb = solve(B[lab], ta * (1 + mp.mpf(10) ** -8))
        print(lab, "B", mp.nstr(tb, 66), " agree %.1f digits" % float(-mp.log10(abs(ta - tb) / ta)) if ta != tb else " agree all", flush=True)
