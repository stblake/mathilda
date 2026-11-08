"""Holdout B, family C4 (experiment 229) -- reference values.

Two independent routes to each global minimum:
  A (this file, mode "A"): scipy differential_evolution, 3 seeds, popsize 40,
    tol 1e-14 -> best point -> mpmath.findroot on the mpmath.diff gradient at
    dps 90 (bc4.6: the minimum is ON the edge x = 512 -- 1-D findroot in y
    there, plus the KKT sign check df/dx < 0).
  B (mode "B"): a different global search -- L-BFGS-B multistart from a
    uniform grid of starts covering the box (sizes below) -- whose best few
    points are written to bc4_starts.gp for PARI/gp (siam100_refs/bc4_ref.gp),
    which polishes them by Newton with derivnum gradients/Hessians at
    realprecision 110.
The two must agree on the minimiser's basin and the value to >= 60 digits.
Run:  python3 bc4_ref.py A   ;   python3 bc4_ref.py B && gp -q bc4_ref.gp
"""
import sys, os, itertools
import numpy as np
import mpmath as mp
from scipy.optimize import differential_evolution, minimize

import math

# the ten objectives, written out afresh here (not imported from the .py half,
# so a slip in a half cannot leak into its truth); lib = math or mpmath
def make(L):
    sin, cos, exp, sqrt = L.sin, L.cos, L.exp, L.sqrt
    pi, e = (math.pi, math.e) if L is math else (mp.pi, mp.e)
    q = (lambda a, b: a / b) if L is math else (lambda a, b: mp.mpf(a) / b)
    def g5(t):
        return sum(i * cos((i + 1) * t + i) for i in range(1, 6))
    def f1(x, y): return sin(10*x) + cos(11*y) + (x**2 + y**2)/4
    def f2(x, y): return (x**2 + y - 11)**2 + (x + y**2 - 7)**2 + 3*x
    def f3(x, y, z):
        return exp(sin(30*x)) + sin(35*exp(y)) + sin(sin(40*z)) - sin(10*(x + y + z)) + (x**2 + y**2 + z**2)/4
    def f4(x, y):
        return -20*exp(-sqrt(((x - q(3, 10))**2 + (y + q(2, 5))**2)/2)/5) - exp((cos(2*pi*x) + cos(2*pi*y))/2) + e + 20
    def f5(x, y): return g5(x)*g5(y)
    def f6(x, y): return -(y + 47)*sin(sqrt(abs(x/2 + y + 47))) - x*sin(sqrt(abs(x - (y + 47))))
    def f7(x): return cos(x) + cos(sqrt(2)*x) + cos(sqrt(3)*x)
    def f8(x, y):
        return (sin(25*x) + sin(25*y))/2 + (x**2 + y**2)/10 - 2*exp(-300*((x - q(3, 5))**2 + (y + q(7, 20))**2))
    def f9(x, y): return 100*(y - x**2)**2 + (1 - x)**2 + q(3, 10)*sin(15*x)*sin(15*y)
    def f10(x, y, z, w):
        return sum(v**2 - 10*cos(2*pi*v) for v in (x, y, z, w)) + (x*y + y*z + z*w)/2 + q(3, 10)*x
    return [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]


FM, FP = make(math), make(mp)
BOX = [[(-1, 1)] * 2, [(-5, 5)] * 2, [(-1, 1)] * 3, [(-5, 5)] * 2, [(-10, 10)] * 2,
       [(-512, 512)] * 2, [(0, 100)], [(-1, 1)] * 2, [(-2, 2)] * 2, [(-5, 5)] * 4]
GRID = [60, 40, 30, 100, 200, 200, 5000, 80, 80, 11]


def route_a(k):
    f, d = FM[k], len(BOX[k])
    best = min((differential_evolution(lambda v: f(*v), BOX[k], tol=1e-14, seed=s, popsize=40)
                for s in range(3)), key=lambda r: r.fun)
    mp.mp.dps = 90
    F = FP[k]
    x0 = [mp.mpf(float(t)) for t in best.x]
    if k == 5:                                       # bc4.6: edge x = 512
        y = mp.findroot(lambda y: mp.diff(lambda t: F(mp.mpf(512), t), y), x0[1])
        v = [mp.mpf(512), y]
        print("  KKT df/dx at edge:", mp.nstr(mp.diff(lambda t: F(t, y), mp.mpf(512)), 5))
    else:
        g = [(lambda i: (lambda *u: mp.diff(F, u, tuple(int(j == i) for j in range(d)))))(i)
             for i in range(d)]
        v = [mp.findroot(g[0], x0[0])] if d == 1 else list(mp.findroot(g, x0))
    val = F(*v)
    print("bc4.%d A de=%.17g at %s\n  value %s" % (k + 1, best.fun, [mp.nstr(t, 20) for t in v],
                                                   mp.nstr(val, 75)), flush=True)


def route_b(k, out):
    f, d, (lo, hi) = FM[k], len(BOX[k]), BOX[k][0]
    g = np.linspace(lo, hi, GRID[k])
    res = []
    for s in itertools.product(g, repeat=d):
        r = minimize(lambda v: f(*v), np.array(s), method="L-BFGS-B", bounds=BOX[k])
        res.append((r.fun, tuple(r.x)))
    res.sort()
    print("bc4.%d B multistart best %.17g at %s (n starts %d)" % (k + 1, res[0][0], res[0][1], len(res)),
          flush=True)
    out.write("starts[%d] = %s;\n" % (k + 1, "[" + ",".join(
        "[" + ",".join("%.17g" % t for t in p) + "]" for _, p in res[:3]) + "]"))


if __name__ == "__main__":
    mode = sys.argv[1]
    ks = [int(a) - 1 for a in sys.argv[2:]] or range(10)
    if mode == "A":
        for k in ks:
            route_a(k)
    else:
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "bc4_starts.gp"), "w") as out:
            out.write("starts = vector(10);\n")
            for k in ks:
                route_b(k, out)
