"""Reference for d7 (area where exp(-(y+x^3)^2) > y^2/32 + exp(sin y)).
For fixed y the set is the x-interval cbrt(-y-s) < x < cbrt(-y+s), s = sqrt(-log g(y)),
wherever g(y) < 1.  Two islands: y in [y1, 0] and [y2, y3].  mpmath tanh-sinh with
breakpoints at the band ends and at yc (where -y-s = 0, a cube-root kink)."""
from mpmath import mp, mpf, exp, sin, log, sqrt, cbrt, findroot, quad, nstr
mp.dps = 90
G = lambda y: y**2/32 + exp(sin(y))
def s(y): return sqrt(-log(G(y)))
def rc(t): return cbrt(t) if t >= 0 else -cbrt(-t)
def w(y):
    g = G(y)
    if g >= 1: return mpf(0)
    q = sqrt(-log(g))
    return rc(-y + q) - rc(-y - q)
y1 = findroot(lambda y: G(y) - 1, mpf(-2.85))
y2 = findroot(lambda y: G(y) - 1, mpf(3.8))
y3 = findroot(lambda y: G(y) - 1, mpf(4.45))
yc = findroot(lambda y: -log(G(y)) - y**2, mpf(-0.834))
print("y1", nstr(y1, 20), "yc", nstr(yc, 20), "y2", nstr(y2, 20), "y3", nstr(y3, 20))
A = quad(w, [y1, yc, 0]) + quad(w, [y2, y3])
print(nstr(A, 80))
