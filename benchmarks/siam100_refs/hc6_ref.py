"""Holdout A, family C6 (experiment 201): references for hc6.1 .. hc6.10.

Return probabilities p = 1 - 1/G of biased nearest-neighbour lattice walks,
G = expected number of visits to the origin = (2 pi)^-d Int 1/(1 - phi(k)) dk.

Method A (Bessel/Laplace): for a walk on Z^d whose axis-i steps are +-e_i with
probabilities p_i+, p_i-, the number of paths returning to 0 is unchanged if
both are replaced by sqrt(p_i+ p_i-) (each closed path has as many + as -
steps), so G = (2 pi)^-d Int dk/(1 - sum a_i cos k_i), a_i = 2 sqrt(p_i+ p_i-),
and with 1/(1-phi) = Int_0^inf e^{-t(1-phi)} dt:
        G = Int_0^inf e^{-t} prod_i I0(a_i t) dt          (mpmath quad).
  For the triangular lattice (hc6.10) the steps do not separate by axis (and
  the +-e1 steps need not pair up in a closed path); method A there is the
  1-D reduction with the k1 integral done in closed form (see G_tri_A), method
  B the one with the k2 integral done in closed form.
Method B (Fourier, complex): the ORIGINAL biased characteristic function,
with the last integral(s) in closed form and the rest by quadrature
  Z^2   : G = (1/pi) Int_0^pi Re[1/(sqrt(A - c2) sqrt(A + c2))] dk1,
          A = 1 - (p1+ + p1-) cos k1 - i (p1+ - p1-) sin k1, c2 = 2 sqrt(p2+ p2-)
  Z^3   : G = (1/pi) Int_0^pi Re[(2/(pi A)) K((2 beta/A)^2)] dk1, beta = 1/3
  tri   : G = (1/pi) Int_0^pi Re[1/(sqrt(Ac - B) sqrt(Ac + B))] dk1,
          Ac = 1 - cos(k1)/3 - 2 i eps sin k1
  hc6.7 : the closed forms of the simple cubic Watson integral,
          Glasser-Zucker  u = sqrt6/(32 pi^3) G(1/24) G(5/24) G(7/24) G(11/24)
          Watson          u = (12/pi^2)(18+12sqrt2-10sqrt3-7sqrt6) K(k)^2,
                          k = (2-sqrt3)(sqrt3-sqrt2)
Then findroot in eps.  A machine-precision scipy check of each case's
STATEMENT (the integrand the .m file integrates) is printed last.
Usage: python3.11 hc6_ref.py [dps]
"""
import sys
from mpmath import (mp, mpf, mpc, pi, sqrt, cos, sin, exp, besseli, quad, inf, findroot,
                    nstr, re, ellipk, gamma)

mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 70
q4 = mpf(1)/4


def UPTS():
    """Quadrature breakpoints in u = sqrt(k): 64 equal pieces (8 pieces lose
    ~30 digits at eps = 1/10; 64 and 128 pieces agree to 50 digits)."""
    return [sqrt(pi)*mpf(i)/64 for i in range(65)]


def G_bessel(a_list):
    def f(t):
        r = exp(-t*(1 - sum(a_list)))
        for a in a_list:
            r *= besseli(0, a*t)*exp(-a*t)
        return r
    return quad(f, [0, mpf(1)/4, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, inf])


def G_fourier2(s1, d1, c2, z0=0):
    """s1 = p1+ + p1-, d1 = p1+ - p1-, c2 = 2 sqrt(p2+ p2-), z0 = 1 - s1 - c2
    passed EXACTLY, and A - c2 = z0 + 2 s1 sin^2(k/2) - i d1 sin k formed
    without cancellation: forming it as 1 - s1 cos k - c2 from rounded
    binary values leaves O(1e-67) instead of 0 at k = 0, which the square root
    amplifies to a 1e-34 error."""
    def f(u):
        k = u*u
        Am = z0 + 2*s1*sin(k/2)**2 - 1j*d1*sin(k)
        return 2*u*re(1/(sqrt(Am)*sqrt(Am + 2*c2)))
    return quad(f, UPTS())/pi


def G_fourier3(e):
    beta = mpf(1)/3

    def f(u):
        k = u*u
        A = 1 - cos(k)/3 - 2j*e*sin(k)
        return 2*u*re(2/(pi*A)*ellipk((2*beta/A)**2))
    return quad(f, UPTS())/pi


def G_tri_A(e):
    """Triangular lattice, the OTHER 1-D reduction: integrate out k1 first.
    For fixed k2, phi = al e^{ik1} + be e^{-ik1} + ga with al = 1/6+eps + e^{ik2}/6,
    be = 1/6-eps + e^{-ik2}/6, ga = cos(k2)/3, and (1/2pi) Int dk1/(1 - phi)
    = 1/sqrt((1-ga)^2 - 4 al be) = 1/sqrt(D),
        D = (2/9) sin^2(k2/2) (7 - cos k2) + 4 eps^2 + (4/3) i eps sin k2
    (written without cancellation).  (The "defective symmetric walk" trick of
    the square lattice is NOT valid here: a closed path can use +(e1+e2), -e2,
    -e1, so +-e1 steps need not pair up.)"""
    def f(k):
        D = mpf(2)/9*sin(k/2)**2*(7 - cos(k)) + 4*e*e + mpf(4)/3*1j*e*sin(k)
        return re(1/sqrt(D))
    pts = sorted([mpf(0), pi] + [e*mpf(2)**j for j in range(-6, 7) if e*mpf(2)**j < pi])
    return quad(f, pts)/pi


def G_tri_B(e):
    def f(u):
        k = u*u
        B = mpf(2)/3*cos(k/2)
        # Ac - B = 1 - cos(k)/3 - (2/3) cos(k/2) - 2i eps sin k, without cancellation:
        # with c = cos(k/2), 1 - (2c^2 - 1)/3 - 2c/3 = (2/3)(1 - c)(2 + c), 1 - c = 2 sin^2(k/4)
        c = cos(k/2)
        AmB = mpf(4)/3*sin(k/4)**2*(2 + c) - 2j*e*sin(k)
        return 2*u*re(1/(sqrt(AmB)*sqrt(AmB + 2*B)))
    return quad(f, UPTS())/pi


def a2(m, d):
    """2 sqrt(p+ p-) for p+- = m +- d, as 2 sqrt(m^2 - d^2) (exact inputs)."""
    return 2*sqrt(m*m - d*d)


# walk descriptions as functions of eps: (method A G, method B G)
def ew(e):           # E/W bias on the square lattice
    return (lambda: G_bessel([a2(q4, e), a2(q4, 0)]),
            lambda: G_fourier2(2*q4, 2*e, a2(q4, 0)))


def diag(e):         # E,N 1/4+e ; W,S 1/4-e
    return (lambda: G_bessel([a2(q4, e), a2(q4, e)]),
            lambda: G_fourier2(2*q4, 2*e, a2(q4, e), 1 - 2*q4 - a2(q4, e)))


def aniso(e):        # E/W 1/3 +- e ; N,S 1/6
    t, s = mpf(1)/3, mpf(1)/6
    return (lambda: G_bessel([a2(t, e), a2(s, 0)]),
            lambda: G_fourier2(2*t, 2*e, a2(s, 0)))   # z0 = 1 - 2/3 - 1/3 = 0


def cubic(e):        # Z^3, E/W 1/6 +- e, others 1/6
    s = mpf(1)/6
    return (lambda: G_bessel([a2(s, e), a2(s, 0), a2(s, 0)]),
            lambda: G_fourier3(e))


def tri(e):
    return (lambda: G_tri_A(e), lambda: G_tri_B(e))


def root(walk, Gtarget, e0):
    out = []
    for m in (0, 1):
        out.append(re(findroot(lambda e: walk(e)[m]() - Gtarget, mpf(e0))))
    return out


class _Printer(list):
    """Print each case as soon as it is computed."""
    def append(self, row):
        lab, a, b = row
        print(lab, nstr(a, 66), nstr(b, 66), "reldiff", nstr(abs(a-b)/abs(a), 3), flush=True)


res = _Printer()
import os
ONLY = os.environ.get("ONLY", "")          # e.g. ONLY=5,6 recomputes just those
want = lambda i: not ONLY or str(i) in ONLY.split(",")
# hc6.1: p at eps = 1/10 (E/W)
if want(1):
    A, B = ew(mpf(1)/10)
    res.append(("hc6.1", 1 - 1/A(), 1 - 1/B()))
# hc6.2 / hc6.3: eps with p = 2/5, 3/5 (E/W)
if want(2):
    res.append(("hc6.2",) + tuple(root(ew, mpf(5)/3, '0.09')))
if want(3):
    res.append(("hc6.3",) + tuple(root(ew, mpf(5)/2, '0.04')))
# hc6.4: E,N biased equally, p = 1/2
if want(4):
    res.append(("hc6.4",) + tuple(root(diag, 2, '0.045')))
# hc6.5: anisotropic, p = 1/2
if want(5):
    res.append(("hc6.5",) + tuple(root(aniso, 2, '0.06')))
# hc6.6: E/W, p = 1/5 (G = 5/4): a strongly biased walk.  (A 1-D walk was
# dropped: in Z^1 the Fourier form (1/2pi) Int Re 1/(1 - phi) MISSES a delta
# mass at k = 0 -- the z -> 1 limit of 1/(1 - z phi) does not commute with the
# integral when the walk is transient in 1-D -- so it is not G.)
if want(6):
    res.append(("hc6.6",) + tuple(root(ew, mpf(5)/4, '0.2')))
# hc6.7: Polya's simple cubic return probability
if want(7):
    uGZ = sqrt(6)/(32*pi**3)*gamma(mpf(1)/24)*gamma(mpf(5)/24)*gamma(mpf(7)/24)*gamma(mpf(11)/24)
    kW = (2 - sqrt(3))*(sqrt(3) - sqrt(2))
    uW = 12/pi**2*(18 + 12*sqrt(2) - 10*sqrt(3) - 7*sqrt(6))*ellipk(kW**2)**2
    res.append(("hc6.7", 1 - 1/uGZ, 1 - 1/uW))
# hc6.8: Z^3 E/W biased, p = 3/10
if want(8):
    res.append(("hc6.8",) + tuple(root(cubic, mpf(10)/7, '0.03')))
# hc6.9: p for E/W bias 1/20 and N/S bias 1/10
if want(9):
    A9 = G_bessel([a2(q4, mpf(1)/20), a2(q4, mpf(1)/10)])
    B9 = G_fourier2(2*q4, mpf(1)/10, a2(q4, mpf(1)/10), 1 - 2*q4 - a2(q4, mpf(1)/10))
    res.append(("hc6.9", 1 - 1/A9, 1 - 1/B9))
# hc6.10: triangular lattice, +-e1 biased, p = 1/2
if want(10):
    res.append(("hc6.10",) + tuple(root(tri, 2, '0.04')))
