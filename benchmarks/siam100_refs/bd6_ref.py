# Holdout B, family D6 (smooth oscillatory towers on a finite interval).
# Method A (this file): mpmath tanh-sinh at dps 80, each interval split at the
# points where the oscillating phase phi(x) passes k*pi (inverse phase given
# in closed form per case).  Method B: bd6_ref.gp (PARI/gp intnum, split at
# the points phi = k*pi + pi/2) -- different software, different split.
# bd6.7 additionally has the closed form pi*BesselJ(3, 50).
import mpmath as mp
mp.mp.dps = 80
E = mp.e
third = mp.mpf(1)/3; half = mp.mpf(1)/2
# (label, f, a, b, phi(x), phi^{-1}(y))
cases = [
 ("bd6.1", lambda x: mp.exp(mp.sin(mp.exp(x))), 0, 1, lambda x: mp.exp(x), lambda y: mp.log(y)),
 ("bd6.2", lambda x: mp.cos(x + mp.exp(mp.exp(x + half))), -1, 1, lambda x: mp.exp(mp.exp(x+half)), lambda y: mp.log(mp.log(y)) - half),
 ("bd6.3", lambda x: mp.exp(mp.cos(x**3 + mp.exp(x))), 0, 2, lambda x: x**3 + mp.exp(x), None),
 ("bd6.4", lambda x: mp.log(2 + mp.sin(mp.exp(mp.exp(x)))), -1, 1, lambda x: mp.exp(mp.exp(x)), lambda y: mp.log(mp.log(y))),
 ("bd6.5", lambda x: 1/(mp.mpf(5)/4 + mp.cos(mp.exp(4*x))), 0, 1, lambda x: mp.exp(4*x), lambda y: mp.log(y)/4),
 ("bd6.6", lambda x: mp.sin(mp.exp(mp.exp(2*x))), 0, 1, lambda x: mp.exp(mp.exp(2*x)), lambda y: mp.log(mp.log(y))/2),
 ("bd6.7", lambda x: mp.cos(50*mp.sin(x) - 3*x), 0, mp.pi, None, None),
 ("bd6.8", lambda x: mp.exp(x)*mp.cos(mp.sinh(mp.exp(x + half))), -1, 1, lambda x: mp.sinh(mp.exp(x+half)), lambda y: mp.log(mp.asinh(y)) - half),
 ("bd6.9", lambda x: mp.cos(mp.exp(mp.exp(mp.exp(x)))), -1, half, lambda x: mp.exp(mp.exp(mp.exp(x))), lambda y: mp.log(mp.log(mp.log(y)))),
 ("bd6.10", lambda x: mp.exp(x + mp.sin(mp.exp(mp.exp(x + half)))), -1, 1, lambda x: mp.exp(mp.exp(x+half)), lambda y: mp.log(mp.log(y)) - half),
]
for lab, f, a, b, phi, inv in cases:
    a = mp.mpf(a); b = mp.mpf(b)
    if phi is None:
        pts = mp.linspace(a, b, 41)
    else:
        lo, hi = phi(a), phi(b)
        ks = range(int(mp.ceil(lo/mp.pi)), int(mp.floor(hi/mp.pi)) + 1)
        if inv is None:
            mids = [mp.findroot(lambda x: phi(x) - k*mp.pi, (a, b), solver='anderson') for k in ks]
        else:
            mids = [inv(k*mp.pi) for k in ks]
        pts = [a] + [m for m in mids if a < m < b] + [b]
    v = mp.quad(f, pts)
    print(lab, mp.nstr(v, 70), len(pts))
print("bd6.7 closed", mp.nstr(mp.pi*mp.besselj(3, 50), 70))
