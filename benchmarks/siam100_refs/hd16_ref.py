"""Reference A for holdout family hd16 (exp 221): heat-equation crossing times.

u_t = u_xx (hd16.10: u_t = u_xx + u_yy), each problem a time t* at which a
functional of u(., t) reaches a level.  Method A (this file): the separated
(Fourier sine/cosine) series of u -- coefficients of the initial data in
closed form -- summed directly to |term| < 1e-80, arc lengths by mpmath.quad
(split at the kinks), the root by mpmath.findroot; hd16.1 uses the closed
form L(t) = (2/pi) sqrt(1+a^2) E(a^2/(1+a^2)), a = pi exp(-pi^2 t).
Reference B (siam100_refs/hd16_ref.wl) uses the METHOD OF IMAGES (erf /
Gaussian image sums, no Fourier coefficients).
Usage: python3.11 hd16_ref.py
"""
import mpmath as mp
mp.mp.dps = int(__import__("os").environ.get("DPS", "75"))
pi = mp.pi
TOL = mp.mpf(10) ** -(mp.mp.dps + 10)


def ssum(term, t, lam):
    """Sum_{k>=1} term(k) e^{-lam(k) t}, to |term| < TOL."""
    if t <= 0:
        raise ValueError("t <= 0")
    s, k = mp.mpf(0), 1
    while True:
        e = mp.exp(-lam(k) * t)
        if e < TOL and k > 3:
            return s
        s += term(k) * e
        k += 1


odd = lambda k: k % 2 == 1
sk2 = lambda k: [0, 1, 0, -1][k % 4]          # sin(k pi/2)
L1 = lambda k: (k * pi) ** 2


def center_hat(t):   return ssum(lambda k: 8 / (pi * k) ** 2 * sk2(k) ** 2, t, L1)
def center_one(t):   return ssum(lambda k: (4 / (k * pi) if odd(k) else 0) * sk2(k), t, L1)
def center_par(t):   return ssum(lambda k: (8 / (k * pi) ** 3 if odd(k) else 0) * sk2(k), t, L1)
def center_ind(t):   return ssum(lambda k: 2 / (k * pi) * (mp.cospi(mp.mpf(k) / 4) - mp.cospi(3 * mp.mpf(k) / 4)) * sk2(k), t, L1)
def end_neumann(t):  return mp.mpf(1) / 2 - ssum(lambda k: 4 / (k * pi) ** 2 if odd(k) else 0, t, L1)
def quarter_inhom(t): return mp.mpf(1) / 4 + ssum(lambda k: 2 * (-1) ** k / (k * pi) * mp.sinpi(mp.mpf(k) / 4), t, L1)


def len_sin(t):
    a = pi * mp.exp(-pi ** 2 * t)
    return 2 / pi * mp.sqrt(1 + a * a) * mp.ellipe(a * a / (1 + a * a))


def len_hat(t):
    ux = lambda x: ssum(lambda k: 8 / (pi * k) ** 2 * sk2(k) * k * pi * mp.cospi(k * x), t, L1)
    return 2 * mp.quad(lambda x: mp.sqrt(1 + ux(x) ** 2), [0, mp.mpf(1) / 4, mp.mpf(1) / 2])


def len_ind3(t):
    lam = lambda k: (k * pi / 3) ** 2
    bk = lambda k: 2 / (k * pi) * (mp.cospi(mp.mpf(k) / 3) - mp.cospi(2 * mp.mpf(k) / 3))
    ux = lambda x: ssum(lambda k: bk(k) * k * pi / 3 * mp.cospi(k * x / 3), t, lam)
    # symmetric about x = 3/2
    return 2 * mp.quad(lambda x: mp.sqrt(1 + ux(x) ** 2), [0, mp.mpf(1) / 2, 1, mp.mpf(5) / 4, mp.mpf(3) / 2])


probs = [
    ("hd16.1", lambda t: len_sin(t) - mp.mpf(3) / 2, mp.mpf("0.05")),
    ("hd16.2", lambda t: center_hat(t) - mp.mpf(1) / 2, mp.mpf("0.05")),
    ("hd16.3", lambda t: center_one(t) - mp.mpf(1) / 2, mp.mpf("0.09")),
    ("hd16.4", lambda t: center_par(t) - mp.mpf(1) / 10, mp.mpf("0.09")),
    ("hd16.5", lambda t: center_ind(t) - mp.mpf(1) / 2, mp.mpf("0.06")),
    ("hd16.6", lambda t: end_neumann(t) - mp.mpf(1) / 3, mp.mpf("0.08")),
    ("hd16.7", lambda t: len_hat(t) - mp.mpf(6) / 5, mp.mpf("0.0991")),
    ("hd16.8", lambda t: len_ind3(t) - 4, mp.mpf("0.01822")),
    ("hd16.9", lambda t: quarter_inhom(t) - mp.mpf(1) / 10, mp.mpf("0.1084")),
    ("hd16.10", lambda t: center_one(t) ** 2 - mp.mpf(1) / 2, mp.mpf("0.06")),
]
import sys
only = sys.argv[1:]
for lab, F, t0 in probs:
    if only and lab not in only:
        continue
    print(lab, mp.nstr(mp.findroot(F, (t0 * mp.mpf("0.97"), t0 * mp.mpf("1.03")), tol=mp.mpf(10) ** -(mp.mp.dps - 5)), 70), flush=True)
