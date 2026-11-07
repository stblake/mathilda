"""Holdout A, family C8 (experiment 203): references for hc8.1 .. hc8.10.

Heat-equation crossing times t*: the time a temperature (a point value or the
mean) of a heat problem with zero initial data and boundary data 1 (or initial
data 1 and zero boundary data) reaches a threshold.

  w(t)    = (4/pi) sum (-1)^k/(2k+1) exp(-(2k+1)^2 pi^2 t/4)      centre of [-1,1], zero ends, initial 1
  W(x,t)  = (4/pi) sum (-1)^k/(2k+1) cos((2k+1) pi x/2) exp(...)   the same at x
  m(t)    = (8/pi^2) sum exp(-(2k+1)^2 pi^2 t/4)/(2k+1)^2          its mean over [-1,1]
  disk    = 1 - sum_n 2/(j_n J1(j_n)) exp(-j_n^2 t)                 centre of the unit disk, boundary 1
  ball    = 1 + 2 sum_{n>=1} (-1)^n exp(-n^2 pi^2 t)                centre of the unit ball, surface 1

Method A: the eigen-series above (mpmath nsum; besseljzero for the disk), findroot.
Method B, independent:
  w, W : the method-of-images (short-time) series
         W(x,t) = sum_{n in Z} (-1)^n (erf((x+1-2n)/(2 sqrt t)) - erf((x-1-2n)/(2 sqrt t)))/2
  m    : numerical inverse Laplace transform (mpmath invertlaplace, Talbot) of
         M(s) = (1 - tanh(sqrt s)/sqrt s)/s
  disk : invertlaplace of U(s) = 1/(s I0(sqrt s))  (the centre value's transform)
  ball : Jacobi's imaginary transformation  u = (2/sqrt(pi t)) sum_{n>=0} exp(-(n+1/2)^2/t)
Usage: python3.11 hc8_ref.py
"""
from mpmath import (mp, mpf, pi, exp, erf, sqrt, cos, nsum, inf, findroot, nstr,
                    besseljzero, besselj, besseli, tanh, invertlaplace)

mp.dps = 80
sg = lambda k: (-1)**int(k)

wA = lambda t: 4/pi*nsum(lambda k: sg(k)/(2*k+1)*exp(-(2*k+1)**2*pi**2*t/4), [0, inf])
WA = lambda x, t: 4/pi*nsum(lambda k: sg(k)/(2*k+1)*cos((2*k+1)*pi*x/2)*exp(-(2*k+1)**2*pi**2*t/4), [0, inf])
mA = lambda t: 8/pi**2*nsum(lambda k: exp(-(2*k+1)**2*pi**2*t/4)/(2*k+1)**2, [0, inf])
diskA = lambda t: 1 - nsum(lambda n: 2/(besseljzero(0, int(n))*besselj(1, besseljzero(0, int(n))))
                           * exp(-besseljzero(0, int(n))**2*t), [1, inf])
ballA = lambda t: 1 + 2*nsum(lambda n: sg(n)*exp(-n**2*pi**2*t), [1, inf])


def WB(x, t):
    s = 2*sqrt(t)
    return nsum(lambda n: sg(n)*(erf((x+1-2*n)/s) - erf((x-1-2*n)/s))/2, [-inf, inf])


wB = lambda t: WB(0, t)
mB = lambda t: invertlaplace(lambda s: (1 - tanh(sqrt(s))/sqrt(s))/s, t, method='talbot')
diskB = lambda t: invertlaplace(lambda s: 1/(s*besseli(0, sqrt(s))), t, method='talbot')
ballB = lambda t: 2/sqrt(pi*t)*nsum(lambda n: exp(-(n+mpf(1)/2)**2/t), [0, inf])

h = mpf(1)/2
cases = [
    # label, function of (w, W, m, disk, ball) implementations -> residual, start
    ("hc8.1",  lambda w, W, m, d, b: lambda t: mpf(5)/4*(1 - w(t)**2) - h, mpf('0.2')),   # square, one side 5, centre 1/2
    ("hc8.2",  lambda w, W, m, d, b: lambda t: 1 - w(t)**2 - h, mpf('0.3')),              # square all sides 1
    ("hc8.3",  lambda w, W, m, d, b: lambda t: w(t) - h, mpf('0.3')),                     # rod midpoint falls to 1/2
    ("hc8.4",  lambda w, W, m, d, b: lambda t: m(t) - h, mpf('0.2')),                     # rod mean falls to 1/2
    ("hc8.5",  lambda w, W, m, d, b: lambda t: 1 - w(t)**3 - h, mpf('0.2')),              # cube all faces 1
    ("hc8.6",  lambda w, W, m, d, b: lambda t: 1 - w(t/4)*w(t) - h, mpf('0.4')),          # 4x2 rectangle all sides 1
    ("hc8.7",  lambda w, W, m, d, b: lambda t: d(t) - h, mpf('0.15')),                    # unit disk
    ("hc8.8",  lambda w, W, m, d, b: lambda t: b(t) - h, mpf('0.1')),                     # unit ball
    ("hc8.9",  lambda w, W, m, d, b: lambda t: 1 - W(h, t)*w(t) - h, mpf('0.2')),         # square, point (1/2,0)
    ("hc8.10", lambda w, W, m, d, b: lambda t: 1 - w(t)**2 - mpf(1)/100, mpf('0.06')),    # square, centre 1/100
]

for lab, mk, t0 in cases:
    ta = findroot(mk(wA, WA, mA, diskA, ballA), t0)
    tb = findroot(mk(wB, WB, mB, diskB, ballB), t0)
    print(lab, nstr(ta, 70), nstr(tb, 70), "reldiff", nstr(abs(ta-tb)/ta, 3), flush=True)
