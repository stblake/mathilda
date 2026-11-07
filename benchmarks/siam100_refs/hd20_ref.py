"""Reference A for holdout family hd20 (exp 225): half-iterates h(h(z)) = f(z).

Each problem: f analytic with f(p) = p, f'(p) = lam (0 < lam != 1); h is the
REGULAR (Schroeder) half-iterate at p, h = Psi(sqrt(lam) sigma(z)), where
sigma(f(z)) = lam sigma(z), sigma(p) = 0, sigma'(p) = 1 and Psi = sigma^-1.
Method A (this file), independent of the power-series statement asked of the
systems:
  closed form  where sigma is elementary (hd20.1-5, hd20.6);
  Koenigs limits otherwise (hd20.7-9):
    repelling lam > 1:  sigma(z) = lim lam^n g^n(z)  (g the local inverse of f),
                        Psi(w)   = lim f^n(w / lam^n)
    attracting lam < 1: sigma(z) = lim lam^-n (f^n(z) - p),
                        Psi(w)   = lim g^n(p + lam^n w)
    (errors O(lam^-n) / O(lam^n); n chosen so the error < 1e-75, and the run
    is repeated at 1.3 n to confirm);
  hd20.10 is a FINITE proxy (Robbins' truncated Abel system for base 2, n=16),
  solved by exact rational Gaussian elimination in Q(log 2) carried as mpf at
  120 digits.
Usage: python3.11 hd20_ref.py
"""
import mpmath as mp
mp.mp.dps = 120
M = mp.mpf

def koenigs_rep(f, g, lam, z0, n):
    s = z0
    for _ in range(n):
        s = g(s)
    s = s * lam ** n                       # sigma(z0)
    w = mp.sqrt(lam) * s
    y = w / lam ** n
    for _ in range(n):
        y = f(y)
    return y

def koenigs_att(f, g, lam, p, z0, n):
    s = z0
    for _ in range(n):
        s = f(s)
    s = (s - p) / lam ** n
    w = mp.sqrt(lam) * s
    y = p + lam ** n * w
    for _ in range(n):
        y = g(y)
    return y

def robbins2(n):
    L = mp.log(2)
    A = [[(k * L) ** j / mp.factorial(j) - (1 if j == k else 0) for k in range(1, n + 1)] for j in range(n)]
    v = mp.lu_solve(mp.matrix(A), mp.matrix([1] + [0] * (n - 1)))
    return mp.findroot(lambda x: -M(1) / 2 + sum(v[k] * x ** (k + 1) for k in range(n)), M('0.5'))

res = {}
res[1] = M(8) / 9                                                 # f = 16z/(4+3z): sigma = 4z/(4-z), h(1/2) = 8/9
res[2] = M(20) / 99                                               # f = tan(4 atan z): h = tan(2 atan z)
res[3] = mp.sin(mp.sqrt(2) * mp.asin(mp.sqrt(M(1) / 10))) ** 2     # f = 4z(1-z)
res[4] = mp.sinh(mp.sqrt(2) * mp.asinh(mp.sqrt(M(1) / 10))) ** 2   # f = 4z(1+z)
res[5] = mp.cos(mp.sqrt(3) * mp.acos(M(9) / 10))                   # T3 at x=1, x0 = 9/10
res[6] = mp.tanh(mp.sqrt(2) * mp.atanh(M(1) / 5))                  # f = 2z/(1+z^2)
for n in (260, 340):
    r7 = koenigs_rep(lambda z: 2 * mp.sin(z), lambda z: mp.asin(z / 2), M(2), M(1) / 5, n)
    r8 = koenigs_rep(lambda z: 2 * z * mp.exp(z), lambda z: mp.lambertw(z / 2).real, M(2), M(1) / 10, n)
    with mp.workdps(300):   # backward iteration amplifies rounding by lam^-n ~ 1e90..1e120
        r9 = +koenigs_att(lambda z: mp.sqrt(2) ** z, lambda y: 2 * mp.log(y, 2), mp.log(2), M(2), M(3) / 2, int(n * 1.6))
    print("n", n, mp.nstr(r7, 80), mp.nstr(r8, 80), mp.nstr(r9, 80), sep="\n  ")
res[7], res[8], res[9] = r7, r8, r9
res[10] = robbins2(16)
for k in range(1, 11):
    print("hd20.%d" % k, mp.nstr(res[k], 75))
