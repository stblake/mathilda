#!/usr/bin/env python3
"""roach_corpus_refs.py -- independent reference values for tests/RoachCorpus.m.

Every record of the Roach corpus carries a numeric reference value.  None of
them comes from Mathilda: each is computed here, with mpmath at 40 digits,
directly from the DEFINITION of the quantity --

  * an integral by adaptive quadrature (mpmath.quad / quadosc), on a contour
    rotated into a half plane of exponential decay when the integrand on the
    real axis only converges conditionally (the rotation is exact by Cauchy's
    theorem; the comment on each such case gives the contour);
  * a pFq special value by mpmath.hyper, which sums the series (with its own
    convergence acceleration at |z| = 1).

The paper's closed form is NOT used for the reference.  Where it is transcribed
it is evaluated as well and printed beside the quadrature value, so that a
transcription error in the corpus (or a typo in the paper) shows up as a
disagreement here rather than as a WRONG in the corpus run.

    python3 tools/roach_corpus_refs.py          # table: id, reference, paper, agree
    python3 tools/roach_corpus_refs.py --emit   # the reference values as Mathilda rules

Sources: K. Roach, "Solving Integrals with the Quantum Computer Algebra
System" (2008): 30 definite integrals and 11 pFq special values; "Hypergeometric
Function Representations" (ISSAC '96) sec. 21 gallery: 12 pFq formulas; "Meijer
G Function Representations" (ISSAC '97) sec. 15 gallery: 12 definite integrals.
The sample parameters are the ones in tests/RoachCorpus.m and must be kept in
step with it.
"""

import sys
import mpmath as mp
from mpmath import mpf, pi, sqrt, exp, log, sin, cos, inf, gamma, besselj, \
    bessely, besseli, besselk, erfc, ei, si, ci, atan, acot, asin, hyper, \
    ellipf, ellipk, ellipe, ellippi, zeta, lerchphi, catalan, struveh, \
    angerj, gammainc, atanh, asinh, mpc

mp.mp.dps = 40

R = {}      # id -> reference value
P = {}      # id -> paper's closed form at the sample point (cross-check only)

Q = lambda f, a, b, **kw: mp.quad(f, [a, b], **kw)
half = mpf(1) / 2

# Paper convention: F(phi_s, k) etc. take sin(amplitude) and the MODULUS k.
# mpmath's ellipf(phi, m) takes the amplitude and the PARAMETER m = k^2.
def F_(s, k): return ellipf(mp.asin(s), k * k)
def E_(s, k): return ellipe(mp.asin(s), k * k)
def Pi_(s, n, k): return ellippi(n, mp.asin(s), k * k)
def K_(k): return ellipk(k * k)
def E0(k): return ellipe(k * k)

# ---------------------------------------------------------------------------
# Roach 2008: the 30 definite integrals
# ---------------------------------------------------------------------------

# R08-01  int_-oo^-3 ((z^2+1)(z^2-4))^(-1/2)/(z-5) dz
R["R08-01"] = Q(lambda z: 1 / (sqrt((z**2 + 1) * (z**2 - 4)) * (z - 5)), -inf, -3)
P["R08-01"] = (F_(sqrt(5) / 3, 1 / sqrt(5)) / mpf(5)**1.5 - K_(1 / sqrt(5)) / mpf(5)**1.5
               + atanh(sqrt(mpf(21) / 26)) / sqrt(546) - atanh(sqrt(mpf(13) / 21)) / sqrt(546)
               + 4 / (21 * mpf(5)**1.5) * Pi_(sqrt(5) / 3, mpf(25) / 21, 1 / sqrt(5))
               - 4 / (21 * mpf(5)**1.5) * ellippi(mpf(25) / 21, mpf(1) / 5))

# R08-02  int_0^oo e^(-al x) J0(be x)^2 dx, al = 3/2, be = 1
al, be = mpf(3) / 2, mpf(1)
R["R08-02"] = Q(lambda x: exp(-al * x) * besselj(0, be * x)**2, 0, inf)
P["R08-02"] = 2 / (pi * sqrt(al**2 + 4 * be**2)) * K_(2 * be / sqrt(al**2 + 4 * be**2))

# R08-03  int_0^a (a^2 - x^2)^(n - 1/2) dx, a = 2, n = 1/3
a, n = mpf(2), mpf(1) / 3
R["R08-03"] = Q(lambda x: (a**2 - x**2)**(n - half), 0, a)
P["R08-03"] = gamma((2 * n + 1) / 2) / gamma(n + 1) * sqrt(pi) * a**(2 * n) / 2

# R08-04  int_0^oo e^(-b x - x^2/(4 a)) dx, a = 2, b = 1/2
a, b = mpf(2), half
R["R08-04"] = Q(lambda x: exp(-b * x - x**2 / (4 * a)), 0, inf)
P["R08-04"] = sqrt(a * pi) * exp(a * b**2) * erfc(b * sqrt(a))

# R08-05  int_0^oo x/(e^(3x) - 1)^(1/3) dx
R["R08-05"] = Q(lambda x: x / mp.cbrt(mp.expm1(3 * x)), 0, inf)
P["R08-05"] = pi * log(3) / mpf(3)**1.5 + pi**2 / 27

# R08-06  int_0^oo x e^(-x n)/(e^x + 1) dx, n = 1/2
n = half
R["R08-06"] = Q(lambda x: x * exp(-x * n) / (exp(x) + 1), 0, inf)
P["R08-06"] = lerchphi(-1, 2, n + 1)

# R08-07  int_0^1 ArcSin[x] x^(2n-1) dx, n = 3/4
n = mpf(3) / 4
R["R08-07"] = Q(lambda x: asin(x) * x**(2 * n - 1), 0, 1)
P["R08-07"] = pi / (4 * n) - sqrt(pi) / (4 * n) * gamma((2 * n + 1) / 2) / gamma(n + 1)

# R08-08  int_0^oo x^(a-1)/Sqrt[x^b + 1] dx, a = 3/4, b = 4
a, b = mpf(3) / 4, mpf(4)
R["R08-08"] = Q(lambda x: x**(a - 1) / sqrt(x**b + 1), 0, inf)
P["R08-08"] = 1 / (b * sqrt(pi)) * gamma((b - 2 * a) / (2 * b)) * gamma(a / b)

# R08-09  int_0^(Pi/4) Cos[x]^(-2a-2) Cos[2x]^a dx, a = 1/3
a = mpf(1) / 3
R["R08-09"] = Q(lambda x: cos(x)**(-2 * a - 2) * cos(2 * x)**a, 0, pi / 4)
P["R08-09"] = sqrt(pi) / 2 * gamma(a + 1) / gamma((2 * a + 3) / 2)

# R08-10  int_0^oo x e^(-2 x n)/(e^x + 1) dx, n = 1/2
n = half
R["R08-10"] = Q(lambda x: x * exp(-2 * x * n) / (exp(x) + 1), 0, inf)
P["R08-10"] = lerchphi(-1, 2, 2 * n + 1)

# R08-11  int_0^oo x e^(-3x)/(e^(-x) + 1) dx
R["R08-11"] = Q(lambda x: x * exp(-3 * x) / (exp(-x) + 1), 0, inf)
P["R08-11"] = pi**2 / 12 - mpf(3) / 4

# R08-12  int_0^oo x^(b-1) (a x + 1)^(-c) dx, a = 2, b = 1/2, c = 3/2
a, b, c = mpf(2), half, mpf(3) / 2
R["R08-12"] = Q(lambda x: x**(b - 1) * (a * x + 1)**(-c), 0, inf)
P["R08-12"] = a**(-b) * gamma(b) * gamma(c - b) / gamma(c)

# R08-13  int_0^oo (Si[a x] - Pi/2) e^(-b x) dx, a = 2, b = 1
a, b = mpf(2), mpf(1)
R["R08-13"] = Q(lambda x: (si(a * x) - pi / 2) * exp(-b * x), 0, inf)
P["R08-13"] = atan(a / b) / b - pi / (2 * b)

# R08-14  int_0^oo x e^(-x)/Sqrt[e^x - 1] dx
R["R08-14"] = Q(lambda x: x * exp(-x) / sqrt(mp.expm1(x)), 0, inf)
P["R08-14"] = pi * log(2) - pi / 2

# R08-15  int_0^oo x^(a-1) (b + x)^(-n) dx, a = 1/2, n = 2, b = 3
a, n, b = half, mpf(2), mpf(3)
R["R08-15"] = Q(lambda x: x**(a - 1) * (b + x)**(-n), 0, inf)
P["R08-15"] = gamma(a) * gamma(n - a) / gamma(n) * b**(a - n)

# R08-16  int_0^oo J1(a x) e^(-x^2 b) dx, a = 1, b = 1/2
a, b = mpf(1), half
R["R08-16"] = Q(lambda x: besselj(1, a * x) * exp(-x**2 * b), 0, inf)
P["R08-16"] = (1 - exp(-a**2 / (4 * b))) / a

# R08-17  int_0^oo e^(-a x) x^n/(b + x) dx, a = 1, b = 2, n = 1/2
a, b, n = mpf(1), mpf(2), half
R["R08-17"] = Q(lambda x: exp(-a * x) * x**n / (b + x), 0, inf)
P["R08-17"] = exp(a * b) * gamma(n + 1) * b**n * gammainc(-n, a * b)

# R08-18  int_0^oo ArcCot[x]/(x + 1) dx
R["R08-18"] = Q(lambda x: acot(x) / (x + 1), 0, inf)
P["R08-18"] = catalan + pi * log(2) / 4

# R08-19  int_0^oo Log[x + 1] Cos[a x]/x dx, a = 1 (oscillatory: quadosc)
a = mpf(1)
R["R08-19"] = mp.quadosc(lambda x: log(x + 1) * cos(a * x) / x, [0, inf], omega=a)
P["R08-19"] = ci(a)**2 / 2 + si(a)**2 / 2 - pi * si(a) / 2 + pi**2 / 8

# R08-20  PV int_0^oo ArcTan[x]/(1 - x^2) dx.  x -> 1/x on (1, oo) folds the
# principal value onto int_0^1 (2 ArcTan[t] - Pi/2)/(1 - t^2) dt, regular at 1.
R["R08-20"] = Q(lambda t: (2 * atan(t) - pi / 2) / (1 - t**2), 0, 1)
P["R08-20"] = -catalan

# R08-21  int_0^oo J_mu(a x)/(b + x) dx, mu = 1/3, a = 1, b = 2.  Rotated onto
# the imaginary axis: J = (H1 + H2)/2; H1 decays in the upper half plane, H2 in
# the lower, and the pole x = -b is on neither quarter-plane, so
#   int_0^oo H1(a x)/(b+x) dx = int_0^(i oo), int_0^oo H2(a x)/(b+x) dx = int_0^(-i oo).
mu, a, b = mpf(1) / 3, mpf(1), mpf(2)
def _h1(z): return besselj(mu, z) + 1j * bessely(mu, z)
def _h2(z): return besselj(mu, z) - 1j * bessely(mu, z)
_up = Q(lambda t: _h1(a * 1j * t) / (b + 1j * t) * 1j, 0, inf)
_dn = Q(lambda t: _h2(-a * 1j * t) / (b - 1j * t) * (-1j), 0, inf)
R["R08-21"] = mp.re((_up + _dn) / 2)
P["R08-21"] = pi * angerj(mu, a * b) / sin(mu * pi) - pi * besselj(mu, a * b) / sin(mu * pi)

# R08-22  int_0^oo x^2 Csch[x] e^(-2 x n) dx, n = 1/4
n = mpf(1) / 4
R["R08-22"] = Q(lambda x: x**2 / mp.sinh(x) * exp(-2 * x * n), 0, inf)
P["R08-22"] = zeta(3, (2 * n + 1) / 2) / 2

# R08-23  int_0^oo x Erfc[a x] e^(-x^2 b) dx, a = 1, b = 1/2
a, b = mpf(1), half
R["R08-23"] = Q(lambda x: x * erfc(a * x) * exp(-x**2 * b), 0, inf)
P["R08-23"] = 1 / (2 * b) - a / (2 * b * sqrt(a**2 + b))

# R08-24  int_0^oo Ei[-x^2] e^(x^2 a^2) dx, a = 1/2
a = half
R["R08-24"] = Q(lambda x: ei(-x**2) * exp(x**2 * a**2), 0, inf)
# The paper prints pi^(3/2)/(2a) - Sqrt[pi] ArcSin[a]/a, which is positive,
# while the integrand is negative (Ei(-x^2) < 0). The quadrature equals
# -Sqrt[pi] ArcSin[a]/a (-> -Sqrt[pi] = int Ei(-x^2) dx as a -> 0); the
# printed form is kept below only to show the disagreement.
P["R08-24"] = pi**1.5 / (2 * a) - sqrt(pi) * asin(a) / a

# R08-25 .. R08-30: the elliptic integrals of page 5
R["R08-25"] = Q(lambda z: 1 / ((z + 14) * sqrt((z - 1) * ((z - 13)**2 + 81))), 2, inf)
P["R08-25"] = (K_(3 / sqrt(10)) / mpf(15)**1.5
               - F_(sqrt(15) / 8, 3 / sqrt(10)) / (2 * mpf(15)**1.5)
               - asin(mpf(3)**1.5 / mpf(2)**3.5) / (5 * sqrt(486)))
R["R08-26"] = Q(lambda z: 1 / ((18 - z) * sqrt((1 - z) * ((16 - z)**2 + 64))), -inf, 0)
P["R08-26"] = (K_(1 / sqrt(17)) / mpf(17)**1.5
               - F_(sqrt(17) / 9, 1 / sqrt(17)) / (2 * mpf(17)**1.5)
               - asin(mpf(1) / 9) / 34)
R["R08-27"] = Q(lambda z: sqrt((z**2 - 1) / (9 - z**2)), 1, 2)
P["R08-27"] = 4 * E_(sqrt(mpf(2) / 3), half) - 2 * F_(sqrt(mpf(2) / 3), half) - sqrt(mpf(5) / 3)
R["R08-28"] = Q(lambda z: 1 / ((z + 5) * sqrt((z - 1) * (5 - z) * ((z - 1)**2 + 9))), 4, 5)
P["R08-28"] = (2 / mpf(15)**1.5 * F_(sqrt(5) / 3, 1 / sqrt(5))
               - asin(mpf(1) / 3) / (10 * mpf(3)**1.5))
R["R08-29"] = Q(lambda z: sqrt(z / ((z - 1) * (z - 2))), 2, 3)
P["R08-29"] = (sqrt(6) + mpf(2)**1.5 * F_(1 / sqrt(2), 1 / sqrt(2))
               - mpf(2)**1.5 * E_(1 / sqrt(2), 1 / sqrt(2)))
R["R08-30"] = Q(lambda z: 1 / ((6 - z) * sqrt((5 - z) * (z - 2) * (z - 3) * z**3)), 3, 4)
P["R08-30"] = (F_(sqrt(3) / 2, sqrt(mpf(2) / 3)) / 12 - E_(sqrt(3) / 2, sqrt(mpf(2) / 3)) / 30
               + Pi_(sqrt(3) / 2, mpf(8) / 9, sqrt(mpf(2) / 3)) / 108 + 1 / (45 * mpf(2)**1.5))

# ---------------------------------------------------------------------------
# Roach 2008: the 11 pFq special values (page 6); z = 3/10 for the z-forms
# ---------------------------------------------------------------------------
fr = lambda p, q: mpf(p) / q
z = fr(3, 10)
R["R08-F01"] = hyper([fr(-26, 7), fr(-26, 7)], [1], -1)
P["R08-F01"] = 95 * gamma(fr(1, 7)) * gamma(fr(5, 14)) / (39 * mpf(2)**fr(2, 7) * pi**1.5) * cos(3 * pi / 14)
R["R08-F02"] = hyper([fr(-15, 4), 1, 1], [fr(-7, 2), 2], 1)
P["R08-F02"] = fr(2252, 665) + 18 * log(8) / 19 + 9 * pi / 19 - 36 * log(2) / 19
R["R08-F03"] = hyper([1, 1, fr(3, 2)], [fr(1, 2), fr(9, 2)], -1)
P["R08-F03"] = fr(1841, 15) - 49 * mpf(2)**1.5 * log(sqrt(2) + 1)
R["R08-F04"] = hyper([fr(-1, 2), fr(-1, 2), fr(5, 2)], [fr(3, 2), fr(7, 2)], -1)
P["R08-F04"] = 805 / (9 * mpf(2)**6.5) - 15 * log(sqrt(2) + 1) / 128
R["R08-F05"] = hyper([fr(-11, 4), fr(-5, 2), fr(-5, 2)], [fr(-7, 2), fr(1, 4)], 1)
P["R08-F05"] = 355 * gamma(fr(1, 4))**2 / (7 * mpf(2)**3.5 * sqrt(pi))
R["R08-F06"] = hyper([half, half], [1], fr(9, 10))
P["R08-F06"] = 2 / pi * K_(3 / sqrt(10))
R["R08-F07"] = hyper([1, 1], [fr(2, 3)], fr(8, 9))
P["R08-F07"] = mpf(3)**1.5 * pi + 9 - 3 * log(3)
R["R08-F08"] = hyper([1, 1], [half], fr(-7, 8))
P["R08-F08"] = fr(8, 15) - 8 * sqrt(7) / mpf(15)**1.5 * asinh(sqrt(7) / mpf(2)**1.5)
R["R08-F09"] = hyper([1, 3], [2, 4, 5], -z)
P["R08-F09"] = (144 * besselj(1, 2 * sqrt(z)) / z**3.5 - 6 * (12 - z**2) / z**3
                - 72 * besselj(0, 2 * sqrt(z)) / z**3)
R["R08-F10"] = hyper([fr(-7, 2)], [fr(-5, 2)], -z)
P["R08-F10"] = ((4 * z**2 - 8 * z**3 - 6 * z + 15) * exp(-z) / 15
                - 8 * sqrt(pi) * z**3.5 * mp.erf(sqrt(z)) / 15)
R["R08-F11"] = hyper([2], [fr(9, 4)], z)
P["R08-F11"] = (5 * (1 - 4 * z) * gammainc(fr(5, 4), z) * exp(z) / (16 * z**1.25)
                - 5 * (1 - 4 * z) * gamma(fr(1, 4)) * exp(z) / (64 * z**1.25) + fr(5, 4))

# ---------------------------------------------------------------------------
# ISSAC '96 gallery: 12 pFq formulas, at z = 3/10
# ---------------------------------------------------------------------------
sz = sqrt(z)
def ber(nu, x): return mp.re(besselj(nu, x * exp(3j * pi / 4)))
def bei(nu, x): return mp.im(besselj(nu, x * exp(3j * pi / 4)))
def erfi(x): return mp.erfi(x)
def L_(nu, x): return mp.struvel(nu, x)
R["I96-01"] = hyper([half], [fr(9, 2)], z)
P["I96-01"] = (-(525 + 280 * z + 140 * z**2) / (128 * z**3) * exp(z)
               + (525 + 630 * z + 420 * z**2 + 280 * z**3) / (256 * z**3.5) * sqrt(pi) * erfi(sz))
R["I96-02"] = hyper([fr(-3, 2)], [fr(-1, 2), half], z)
P["I96-02"] = (1 + 2 * z) * mp.cosh(2 * sz) + sz * mp.sinh(2 * sz) - 4 * z**1.5 * mp.shi(2 * sz)
R["I96-03"] = hyper([fr(-3, 2), fr(-1, 2)], [fr(-5, 2), 1], z)
P["I96-03"] = (5 - 4 * z) / 5 * exp(z / 2) * besseli(0, z / 2) + 4 * z / 5 * exp(z / 2) * besseli(1, z / 2)
R["I96-04"] = hyper([fr(-3, 2), fr(-1, 2)], [2], z)
P["I96-04"] = (-(4 + 24 * z - 28 * z**2) / (15 * z * pi) * K_(sz)
               + (4 + 56 * z + 4 * z**2) / (15 * z * pi) * E0(sz))
R["I96-05"] = hyper([fr(-1, 2), 1], [fr(1, 4), half, fr(3, 4)], z)
_z4 = z**fr(1, 4)
P["I96-05"] = (1 + _z4 * sqrt(2) * sqrt(pi) * exp(2 * sz) * mp.erf(sqrt(2) * _z4)
               - _z4 * sqrt(2) * sqrt(pi) * exp(-2 * sz) * erfi(sqrt(2) * _z4)
               - 2 * sz * pi * mp.erf(sqrt(2) * _z4) * erfi(sqrt(2) * _z4))
R["I96-06"] = hyper([], [1, fr(5, 2), 4], z)
_w = 2 * sqrt(2) * _z4
P["I96-06"] = (-18 / z**2 * ber(0, _w)**2
               - (36 - 36 * sz + 9 * z) / (2 * z**fr(9, 4)) * ber(0, _w) * ber(1, _w)
               + (36 + 36 * sz + 9 * z) / (2 * z**fr(9, 4)) * ber(0, _w) * bei(1, _w)
               - (18 + 27 * z) / (2 * z**2.5) * ber(1, _w)**2
               - (36 + 36 * sz + 9 * z) / (2 * z**fr(9, 4)) * ber(1, _w) * bei(0, _w)
               - 18 / z**2 * bei(0, _w)**2
               - (36 - 36 * sz + 9 * z) / (2 * z**fr(9, 4)) * bei(0, _w) * bei(1, _w)
               - (18 + 27 * z) / (2 * z**2.5) * bei(1, _w)**2)
R["I96-07"] = hyper([fr(3, 2)], [fr(5, 2), 5], z)
P["I96-07"] = (-(432 - 24 * z + 96 * z**2) / (5 * z**3) * besseli(0, 2 * sz)
               + (432 + 192 * z + 48 * z**2) / (5 * z**3.5) * besseli(1, 2 * sz)
               - 48 / (5 * z) * pi * (besseli(0, 2 * sz) * L_(1, 2 * sz) - besseli(1, 2 * sz) * L_(0, 2 * sz)))
R["I96-08"] = hyper([fr(-1, 2), 1, 2], [3, 4], z)
P["I96-08"] = (-(480 + 3472 * z - 2100 * z**2) / (525 * z**3)
               + (480 + 3712 * z - 1024 * z**2 + 192 * z**3) / (525 * z**3) * sqrt(1 - z)
               - 32 / (5 * z**2) * log(half + sqrt(1 - z) / 2))
R["I96-09"] = hyper([fr(-1, 2), half, 1], [fr(3, 2), fr(5, 2)], z)
P["I96-09"] = (-(3 - 3 * z) / (16 * z) - (3 - 3 * z**2) / (32 * z**1.5) * (log(1 - sz) - log(1 + sz))
               + 3 / (8 * sz) * mp.polylog(2, sz) - 3 / (8 * sz) * mp.polylog(2, -sz))
R["I96-10"] = hyper([1, 2, 3], [half, 4], z)
P["I96-10"] = (-(45 - 30 * z - 3 * z**2) / (4 * z**2 * (1 - z))
               - (45 - 60 * z + 9 * z**2) / (2 * z**2.5 * (1 - z)**1.5) * asin(sz)
               + 45 / (4 * z**3) * asin(sz)**2)
R["I96-11"] = hyper([fr(1, 3)], [fr(-2, 3), fr(-1, 2), half, 1], z)
P["I96-11"] = ((1 + 3 * sz) / 2 * besseli(0, 4 * _z4) + (1 - 3 * sz) / 2 * besselj(0, 4 * _z4)
               - 7 * _z4 / 4 * besseli(1, 4 * _z4) + 7 * _z4 / 4 * besselj(1, 4 * _z4))
R["I96-12"] = hyper([half], [fr(-1, 3), fr(-1, 2), fr(-2, 3)], z)
_z3 = z**fr(1, 3)
P["I96-12"] = ((2 - 6 * _z3 + 9 * _z3**2 - 18 * z) / 6 * exp(3 * _z3)
               + (4 + 6 * _z3 - 9 * _z3**2 - 36 * z) / 6 * exp(-3 * _z3 / 2) * cos(3 * sqrt(3) * _z3 / 2)
               + _z3 * (2 + 3 * _z3) / 2 * sqrt(3) * exp(-3 * _z3 / 2) * sin(3 * sqrt(3) * _z3 / 2))

# ---------------------------------------------------------------------------
# ISSAC '97 gallery: 12 definite integrals
# ---------------------------------------------------------------------------

# I97-01  int_0^oo Sin[b x]/(x^2 + c^2)^rho dx, b = 1, c = 2, rho = 3/4
b, c, rho = mpf(1), mpf(2), mpf(3) / 4
R["I97-01"] = mp.quadosc(lambda x: sin(b * x) / (x**2 + c**2)**rho, [0, inf], omega=b)

# I97-02  int_0^oo Cos[a x^2 + 2 b x] dx, a = 1, b = 1/2.  Re of
# int_0^oo e^(i(a x^2 + 2 b x)) dx on the ray x = t e^(i Pi/4), where the
# quadratic phase becomes the Gaussian e^(-a t^2) (the arc at infinity
# vanishes in the sector 0 < arg x < Pi/4).
a, b = mpf(1), half
_w = exp(1j * pi / 4)
R["I97-02"] = mp.re(Q(lambda t: exp(1j * (a * (_w * t)**2 + 2 * b * _w * t)) * _w, 0, inf))
P["I97-02"] = (sqrt(pi) * sqrt(2) * cos(b**2 / a) / (4 * sqrt(a))
               + sqrt(pi) * sqrt(2) * sin(b**2 / a) / (4 * sqrt(a))
               - sqrt(pi) * sqrt(2) * cos(b**2 / a) / (2 * sqrt(a)) * mp.fresnelc(sqrt(2) * b / (sqrt(pi) * sqrt(a)))
               - sqrt(pi) * sqrt(2) * sin(b**2 / a) / (2 * sqrt(a)) * mp.fresnels(sqrt(2) * b / (sqrt(pi) * sqrt(a))))

# I97-03  int_0^oo Cos[b x] ArcTan[a/x^2] dx, a = 1, b = 1
a, b = mpf(1), mpf(1)
R["I97-03"] = mp.quadosc(lambda x: cos(b * x) * atan(a / x**2), [0, inf], omega=b)
P["I97-03"] = pi * sin(sqrt(2) * b * sqrt(a) / 2) / b * exp(-sqrt(2) * b * sqrt(a) / 2)

# I97-04  int_0^oo Cos[b Sqrt[x]] J0(c x)/Sqrt[x] dx, b = 1, c = 1.  x = t^2:
# 2 int_0^oo Cos[b t] J0(c t^2) dt.  Rotating t onto a ray does not decay for
# both factors, so use instead the absolutely convergent Mellin-Parseval
# representation: J0(c t^2) = (1/2pi i) int Gamma(s)/Gamma(1-s) (c t^2/2)^(-2s) ds
# is awkward too; the integrand is summed directly with quadosc on the
# zeros of the slower-decaying factor, then cross-checked against the paper.
b, c = mpf(1), mpf(1)
R["I97-04"] = 2 * mp.quadosc(lambda t: cos(b * t) * besselj(0, c * t * t), [0, inf],
                             zeros=lambda k: sqrt(k * pi / c))
_y = b**2 / (8 * c)
P["I97-04"] = (b * pi / (4 * c) * besselj(fr(7, 4), _y)**2 + 36 * pi * c / b**3 * besselj(fr(3, 4), _y)**2
               - b * pi / (4 * c) * besselj(fr(1, 4), _y)**2
               - 6 * pi / b * besselj(fr(3, 4), _y) * besselj(fr(7, 4), _y))

# I97-05  int_0^oo J1(b x) J1(c x)/x^2 dx, b = 1, c = 2
b, c = mpf(1), mpf(2)
R["I97-05"] = mp.quadosc(lambda x: besselj(1, b * x) * besselj(1, c * x) / x**2, [0, inf], omega=c - b)
P["I97-05"] = (2 * b**2 - 2 * c**2) / (3 * b * pi) * K_(b / c) + (2 * b**2 + 2 * c**2) / (3 * b * pi) * E0(b / c)

# I97-06  int_0^oo x^(3/2) J0(b x) J0(c x) dx, b = 2, c = 1.  Divergent as a
# Lebesgue integral (amplitude x^(1/2)); its Abel / Mellin value is the
# analytic continuation of Weber-Schafheitlin in lambda (int J0 J0 x^-lambda,
# lambda = -3/2, smaller argument c):
#   Gamma((1-lambda)/2) c^0 / (2^lambda b^(1-lambda) Gamma((1+lambda)/2))
#     2F1((1-lambda)/2, (1-lambda)/2; 1; c^2/b^2)
b, c, lam = mpf(2), mpf(1), mpf(-3) / 2
R["I97-06"] = (gamma((1 - lam) / 2) / (2**lam * b**(1 - lam) * gamma((1 + lam) / 2))
               * hyper([(1 - lam) / 2, (1 - lam) / 2], [1], c**2 / b**2))
_kk = sqrt(2 - 2 * sqrt((b**2 - c**2) / b**2)) / sqrt(2)
P["I97-06"] = (1 / (2 * sqrt(b) * gamma(fr(3, 4))**2 * (b**2 - c**2)) * K_(_kk)
               + 1 / (b**2.5 * gamma(fr(3, 4))**2 * ((b**2 - c**2) / b**2)**1.5) * K_(_kk)
               - 2 / (b**2.5 * gamma(fr(3, 4))**2 * ((b**2 - c**2) / b**2)**1.5) * E0(_kk))

# I97-07  int_0^oo e^(i p x) H0^(1)(c x) dx, p = 1/2, c = 1.  Both factors decay
# in the upper half plane; on x = i t, H0^(1)(i c t) = (2/(i pi)) K0(c t):
#   = (2/pi) int_0^oo e^(-p t) K0(c t) dt.
p, c = half, mpf(1)
R["I97-07"] = 2 / pi * Q(lambda t: exp(-p * t) * besselk(0, c * t), 0, inf)
P["I97-07"] = 1 / (c * sqrt((c**2 - p**2) / c**2)) - 2 / (pi * c * sqrt((c**2 - p**2) / c**2)) * asin(p / c)

# I97-08  int_0^oo e^(-p x) I1(c x)^2 dx, c = 1, p = 3
c, p = mpf(1), mpf(3)
R["I97-08"] = Q(lambda x: exp(-p * x) * besseli(1, c * x)**2, 0, inf)
P["I97-08"] = -p / (pi * c**2) * E0(2 * c / p) - (2 * c**2 - p**2) / (p * pi * c**2) * K_(2 * c / p)

# I97-09  int_0^oo e^(-p x) I1(c x)^2/x dx, c = 1, p = 3
R["I97-09"] = Q(lambda x: exp(-p * x) * besseli(1, c * x)**2 / x, 0, inf)
P["I97-09"] = -half + p**2 / (2 * pi * c**2) * E0(2 * c / p) + (4 * c**2 - p**2) / (2 * pi * c**2) * K_(2 * c / p)

# I97-10  int_0^oo x^(nu+1) Sin[c x^2/(2a)] K_nu(c x) dx, nu = 1/3, a = 1, c = 1
nu, a, c = mpf(1) / 3, mpf(1), mpf(1)
# sin(x^2) oscillates ever faster: one tanh-sinh pass over (0, oo) loses
# digits (2e-7 here), so the interval is cut every 1/2 up to 80, where the
# K_nu(c x) factor is below 10^-34.
_cuts = mp.linspace(0, 80, 161) + [inf]
R["I97-10"] = mp.quad(lambda x: x**(nu + 1) * sin(c * x**2 / (2 * a)) * besselk(nu, c * x), _cuts)

# I97-11  int_0^oo Sin[b x^2] K_nu(c x) dx, nu = 1/3, b = 1, c = 1
nu, b, c = mpf(1) / 3, mpf(1), mpf(1)
R["I97-11"] = mp.quad(lambda x: sin(b * x**2) * besselk(nu, c * x), _cuts)
_y = c**2 / (8 * b)
P["I97-11"] = (-pi**1.5 / sin(pi * (nu + 1) / 4) * sin(_y) / sin(pi * nu / 2) / (16 * sqrt(b)) * besselj(nu / 2, _y)
               - pi**1.5 / cos(pi * (nu + 1) / 4) * cos(_y) / sin(pi * nu / 2) / (16 * sqrt(b)) * besselj(nu / 2, _y)
               + pi**1.5 / sin(-pi * (nu - 1) / 4) * sin(_y) / sin(pi * nu / 2) / (16 * sqrt(b)) * besselj(-nu / 2, _y)
               + pi**1.5 / cos(-pi * (nu - 1) / 4) * cos(_y) / sin(pi * nu / 2) / (16 * sqrt(b)) * besselj(-nu / 2, _y))

# I97-12  int_0^oo x J_-nu(b x) (Y_nu(c x) - H_nu(c x)) dx, nu = 1/4, b = 1, c = 2.
# mpmath's StruveH at 40 digits is unusable at the large arguments an
# oscillatory quadrature needs, so Y - H goes through its Laplace form
#   H_nu(z) - Y_nu(z) = 2 (z/2)^nu/(Sqrt[pi] Gamma[nu + 1/2])
#                         int_0^oo e^(-z t) (1 + t^2)^(nu - 1/2) dt    (Re z > 0)
# and the x integral, done first, is the Laplace transform of x^(1+nu) J_-nu(b x):
#   int_0^oo e^(-p x) x^(1+nu) J_-nu(b x) dx
#     = (b/2)^-nu/(p^2 Gamma[1 - nu]) 2F1(1, 3/2; 1 - nu; -b^2/p^2),  p = c t.
# The t integrand is finite at 0 and O(t^(2 nu - 3)) at infinity.
nu, b, c = mpf(1) / 4, mpf(1), mpf(2)
_K = 2 * (c / 2)**nu / (sqrt(pi) * gamma(nu + half))
def _lap(t):
    p = c * t
    return (b / 2)**(-nu) / (p**2 * gamma(1 - nu)) * hyper([1, mpf(3) / 2], [1 - nu], -b**2 / p**2)
R["I97-12"] = -_K * Q(lambda t: (1 + t**2)**(nu - half) * _lap(t), 0, inf)
P["I97-12"] = (-2 * b**(-nu - 2) * c**nu * b**2 * cos(pi * nu) / (pi * (b**2 - c**2))
               + 2 * b**(-3 - nu) * c**(nu + 1) * b**2 * cos(pi * nu) / (pi * (b**2 - c**2)))


def main():
    emit = "--emit" in sys.argv
    order = sorted(R, key=lambda k: (k[:3], k))
    if emit:
        for k in order:
            v = R[k]
            if isinstance(v, mpc) and abs(mp.im(v)) > mpf(10)**-30:
                print('  "%s" -> %s + %s I,' % (k, mp.nstr(mp.re(v), 30), mp.nstr(mp.im(v), 30)))
            else:
                print('  "%s" -> %s,' % (k, mp.nstr(mp.re(v), 30)))
        return
    bad = 0
    for k in order:
        v = R[k]
        pv = P.get(k)
        if pv is None:
            print("%-8s %s" % (k, mp.nstr(v, 25)))
            continue
        err = abs(v - pv) / (1 + abs(v))
        ok = err < mpf(10)**-15
        bad += not ok
        print("%-8s %s  paper %s  %s" % (k, mp.nstr(v, 25), mp.nstr(pv, 25), "ok" if ok else "DIFFERS %s" % mp.nstr(err, 3)))
    print("%d references, %d paper forms differ" % (len(R), bad))
    print("expected to differ (see tests/RoachCorpus.m): R08-01 (real part agrees: the\n"
          "paper's Pi with n > 1 is the principal value), R08-24, R08-30, I96-10, I97-06")


if __name__ == "__main__":
    main()
