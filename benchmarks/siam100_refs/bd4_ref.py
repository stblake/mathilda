"""Holdout B, family D4 reference, method A (mpmath, dps 90): roots of
functions given by special functions or by (slowly convergent / cancelling)
series.  Method B (independent formulas, PARI/gp): bd4_ref.gp; for bd4.4 and
bd4.9 method B is the zeta-zero (Mellin-Barnes residue) expansion of R(e^-u)
below, which never touches the Gram / Moebius series.

  bd4.1  li(x) = 0                              (Soldner's constant)
  bd4.2  zeta(s) = 2, s real                    (mpmath.zeta)
  bd4.3  sum_{n>=2} n^-s / log n = 1            (the series itself: head + Euler-Maclaurin)
  bd4.4  log x with R(x) = 10^-3                (mpmath.riemannr)
  bd4.5  log x with R(x) = 1/10                 (mpmath.riemannr)
  bd4.6  smallest x > 0 with sum_{k>=1} mu(k) x^k = 0  (direct sum, k <= 4000)
  bd4.7  first positive zero of the Riemann-Siegel Z function (mpmath.siegelz)
  bd4.8  P(s) = 1, P the prime zeta function    (mpmath.primezeta)
  bd4.9  log x with R(x) = 10^-9                (mpmath.riemannr)
  bd4.10 zeta'(s) = -1, s real > 1              (mpmath.zeta(s, derivative=1))
"""
import sys
import mpmath as mp
from sympy import mobius

mp.mp.dps = 90
out = {}
out["bd4.1"] = mp.findroot(mp.li, 1.45)
out["bd4.2"] = mp.findroot(lambda s: mp.zeta(s) - 2, 1.7)


def S3(s, N=2000):
    f = lambda x: x ** -s / mp.log(x)
    S = mp.fsum(f(mp.mpf(n)) for n in range(2, N)) + f(mp.mpf(N)) / 2
    S += mp.e1((s - 1) * mp.log(N))                  # exact: int_N^oo x^-s/log x dx = E1((s-1) log N)
    c = mp.taylor(f, mp.mpf(N), 41)
    for k in range(1, 21):
        S -= mp.bernoulli(2 * k) / (2 * k) * c[2 * k - 1]
    return S


out["bd4.3"] = mp.findroot(lambda s: S3(s) - 1, 1.6)
out["bd4.4"] = mp.findroot(lambda t: mp.riemannr(mp.exp(t)) - mp.mpf(10) ** -3, -28)
out["bd4.5"] = mp.findroot(lambda t: mp.riemannr(mp.exp(t)) - mp.mpf(1) / 10, -6)
MU = [0] + [int(mobius(k)) for k in range(1, 4001)]
f6 = lambda x: mp.fsum(MU[k] * x ** k for k in range(1, 4001))
out["bd4.6"] = mp.findroot(f6, 0.6)
out["bd4.7"] = mp.findroot(mp.siegelz, 14)
out["bd4.8"] = mp.findroot(lambda s: mp.primezeta(s) - 1, 1.4)
if "skip9" not in sys.argv:
    mp.mp.dps = 90
    out["bd4.9"] = mp.findroot(lambda t: mp.riemannr(mp.exp(t)) - mp.mpf(10) ** -9, -2800)
out["bd4.10"] = mp.findroot(lambda s: mp.zeta(s, derivative=1) + 1, 1.9)
for k, v in out.items():
    print(k, mp.nstr(v, 80), flush=True)
