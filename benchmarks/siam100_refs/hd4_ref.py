"""Holdout A, family hd4 (experiment 209): roots of number-theoretic special
functions, method B (mpmath, dps 90).  Method A: hd4_ref.wl (wolframscript,
WL's own RiemannR / LogIntegral / Zeta / ZetaZero at WP 90).
  1  li(x) = 0                          (Soldner's constant)
  2  R(x) = 100                          R = mpmath.riemannr
  3  R(x) = 2
  4  R(x) = 1000
  5  R(x) = 1/2, 0 < x < 1
  6  sum mu(k) k^-s = 1/zeta(s) = 1/2   -> zeta(s) = 2
  7  sum Lambda(n) n^-s = -zeta'(s)/zeta(s) = 1
  8  first zero of the Riemann-Siegel Z function (= Im of the first zeta zero)
  9  R(x) = 10^10
  10 the SECOND largest positive zero of R: the zeta-zero (Mellin-Barnes)
     expansion of R(e^-u) of siam100_refs/d4_ref.py, scanned for sign changes
     beyond u* = 34142.128 (the largest zero), refined by findroot."""
import mpmath as mp
mp.mp.dps = 90
out = {}
out[1] = mp.findroot(mp.li, 1.45)
out[2] = mp.findroot(lambda x: mp.riemannr(x) - 100, 540)
out[3] = mp.findroot(lambda x: mp.riemannr(x) - 2, 3)
out[4] = mp.findroot(lambda x: mp.riemannr(x) - 1000, 7900)
out[5] = mp.findroot(lambda x: mp.riemannr(x) - mp.mpf(1) / 2, 0.3)
out[6] = mp.findroot(lambda s: mp.zeta(s) - 2, 1.7)
out[7] = mp.findroot(lambda s: -mp.zeta(s, derivative=1) / mp.zeta(s) - 1, 1.7)
out[8] = mp.findroot(mp.siegelz, 14.1)
out[9] = mp.findroot(lambda x: mp.riemannr(x) - mp.mpf(10) ** 10, 2.3e11)
print("check zetazero(1):", mp.nstr(mp.zetazero(1).imag - out[8], 5))
for k in sorted(out):
    print(k, mp.nstr(out[k], 75), flush=True)
# 10: second largest zero, method of d4_ref.py
NZ = 70
zz = [(mp.zetazero(j), None) for j in range(1, NZ + 1)]
zz = [(r, mp.zeta(r, derivative=1)) for r, _ in zz]


def Ru(u):
    s = mp.mpf(0); n = 1
    while True:
        term = (-1) ** (n + 1) * 2 * (2 * mp.pi) ** (2 * n) / ((2 * n + 1) * mp.zeta(2 * n + 1) * u ** (2 * n + 1))
        s += term
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 20):
            break
        n += 1
    for rho, zp in zz:
        s -= 2 * mp.re(mp.gamma(1 - rho) * u ** (rho - 1) / ((1 - rho) * zp))
    return s


prev = None; roots = []
u = mp.mpf(34000)
while u < 80000 and len(roots) < 3:
    v = Ru(u)
    if prev is not None and (v > 0) != (prev[1] > 0):
        roots.append(mp.findroot(Ru, (prev[0], u), solver='anderson'))
        print("sign change in", float(prev[0]), float(u), "root u =", mp.nstr(roots[-1], 30), flush=True)
    prev = (u, v); u += 50
u2 = roots[1]
u2 = mp.findroot(Ru, u2, tol=mp.mpf(10) ** -150)
print("10 t2 =", mp.nstr(-u2, 80))
print("10", mp.nstr(mp.exp(-u2), 75))
