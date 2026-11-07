"""Reference for holdout-B family d14 (experiment 249): a finite-element
boundary-value problem inside a parameter search (find the parameter for which
a functional of the PDE solution takes a prescribed value).

Every problem is chosen to have a solution in closed form or a rapidly
convergent series, so the truth is known to >= 55 digits.  Two independent
methods per problem (A and B below; each line says which).
Usage: bd14_ref.py k [dps]
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import (mp, mpf, exp, log, sqrt, ei, euler, quad, findroot, nstr, pi, sech, tanh, besseli,
                    besseljzero, nsum, inf, erfi, odefun, dirichlet, cosh)
from bhold_ode_util import bs_solve

k = int(sys.argv[1]); mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 70
out = {}
t0 = time.time()
if k == 1:   # unit disk, -Lap u = exp(a r^2), u = 0 on r = 1: u(0) = 1/2
    A = findroot(lambda a: (ei(a) - euler - log(a)) / (4 * a) - mpf(1) / 2, mpf(2.5))
    B = findroot(lambda a: quad(lambda r: -r * log(r) * exp(a * r * r), [0, 1]) - mpf(1) / 2, mpf(2.5))
    out = {"A u(0) = (Ei(a) - gamma - log a)/(4a)": A, "B Green's function u(0) = Int_0^1 -r log r e^(a r^2) dr": B}
elif k == 2:  # [0,L]x[0,1], -Lap u = 1, u = 0 on x=0, x=L, y=0, du/dy = 0 on y=1: Int_0^L u(x,1) dx = 1/2
    def IA(L):   # sine series in x, cosh in y; the sech-free part summed in closed form (L^3/12)
        return L ** 3 / 12 - nsum(lambda j: 8 * L ** 3 / (pi ** 4 * (2 * j + 1) ** 4) * sech((2 * j + 1) * pi / L), [0, inf])
    def IB(L):   # sine series in y (mu = (m-1/2) pi), cosh in x; closed part L/2 - (64/pi^4) beta(4)
        beta4 = dirichlet(4, [0, 1, 0, -1])
        return L / 2 - 64 * beta4 / pi ** 4 + nsum(lambda m: (-1) ** (m + 1) * 4 / ((m - mpf(1) / 2) * pi) ** 4
                                               * (1 - tanh((m - mpf(1) / 2) * pi * L / 2)), [1, inf])
    A = findroot(lambda L: IA(L) - mpf(1) / 2, mpf(2.3))
    B = findroot(lambda L: IB(L) - mpf(1) / 2, mpf(2.3))
    out = {"A x-sine series (L^3/12 - sum 8L^3 sech(n pi/L)/(pi n)^4)": A,
           "B y-eigen series (L/2 - 64 beta(4)/pi^4 + sum ... (1 - tanh))": B}
elif k == 3:  # equilateral triangle of side s, -Lap u = 1, u = 0: Int u = 1/100
    A = (mpf(320) / 100 / sqrt(3)) ** (mpf(1) / 4)
    def Iq(s):
        u = lambda x, y: y * (sqrt(3) * x - y) * (sqrt(3) * (s - x) - y) / (2 * sqrt(3) * s)
        return 2 * quad(lambda x: quad(lambda y: u(x, y), [0, sqrt(3) * x]), [0, s / 2])
    B = findroot(lambda s: Iq(s) - mpf(1) / 100, mpf(1.17))
    out = {"A Int u = sqrt(3) s^4/320 (exact cubic torsion function)": A,
           "B 2-D quadrature of u = y(sqrt3 x - y)(sqrt3(s-x) - y)/(2 sqrt3 s)": B}
elif k == 4:  # annulus 1/2<r<1, -Lap u = exp(a r^2), u(1) = 0, du/dr(1/2) = 0: Int_{r=1/2} u ds = e^a
    def uA(a):
        return (ei(a) - ei(a / 4)) / (4 * a) - exp(a / 4) * log(2) / (2 * a)
    A = findroot(lambda a: pi * uA(a) - exp(a), mpf(-2))
    def uB(a):   # radial ODE from r = 1/2 (regular point): w(1/2) = 0, w'(1/2) = 0, u(1/2) = -w(1)
        f = odefun(lambda r, w: [w[1], -w[1] / r - exp(a * r * r)], mpf(1) / 2, [mpf(0), mpf(0)])
        return -f(1)[0]
    B = findroot(lambda a: pi * uB(a) - exp(a), A + mpf(10) ** -8)
    out = {"A u(1/2) = (Ei(a) - Ei(a/4))/(4a) - e^(a/4) log 2/(2a)": A, "B radial ODE (odefun) from the Neumann circle": B}
elif k == 5:  # unit disk Bratu -Lap u = lam e^u, u = 0: u(0) = 1 (lower branch)
    A = 8 * (exp(mpf(1) / 2) - 1) / exp(1)
    r0 = mpf(10) ** -12
    def u1(lam):  # radial shooting from r0 with the series u = 1 - (lam e/4) r^2 + (lam e)^2 r^4/64 + O(r^6)
        c = lam * exp(1)
        f = odefun(lambda r, w: [w[1], -w[1] / r - lam * exp(w[0])], r0, [1 - c * r0 ** 2 / 4 + c * c * r0 ** 4 / 64, -c * r0 / 2 + c * c * r0 ** 3 / 16])
        return f(1)[0]
    B = findroot(u1, mpf(1.9))
    out = {"A Liouville: lam = 8 (e^(1/2) - 1)/e": A, "B radial odefun shooting": B}
elif k == 6:  # unit disk -Lap u + kappa^2 u = 1, u = 0: Int u = pi/10
    A = findroot(lambda q: (1 - 2 * besseli(1, q) / (q * besseli(0, q))) / q ** 2 - mpf(1) / 10, mpf(1.5))
    def IB(q):
        return nsum(lambda n: 4 / (besseljzero(0, int(n)) ** 2 * (besseljzero(0, int(n)) ** 2 + q * q)), [1, inf])
    B = findroot(lambda q: IB(q) - mpf(1) / 10, mpf(1.5))
    out = {"A Int u = pi (1 - 2 I1/(k I0))/k^2": A, "B Fourier-Bessel series sum 4 pi/(j^2 (j^2 + k^2))": B}
elif k == 7:  # [0,a]x[0,1] Dirichlet, -Lap u = 1: u(centre) = 1/10
    sA = lambda a: a * a / 8 - nsum(lambda j: (-1) ** j * 4 * a * a / (pi * (2 * j + 1)) ** 3 * sech((2 * j + 1) * pi / (2 * a)), [0, inf])
    sB = lambda a: mpf(1) / 8 - nsum(lambda j: (-1) ** j * 4 / (pi * (2 * j + 1)) ** 3 * sech((2 * j + 1) * pi * a / 2), [0, inf])
    A = findroot(lambda a: sA(a) - mpf(1) / 10, mpf(1.4))
    B = findroot(lambda a: sB(a) - mpf(1) / 10, mpf(1.4))
    out = {"A x-sine series a^2/8 - ...sech(n pi/(2a))": A, "B y-sine series 1/8 - ...sech(n pi a/2)": B}
elif k == 8:  # unit ball -Lap u = exp(a r^2), u = 0: u(0) = 1/4
    def uA(a):
        return (exp(a) - 1) / (2 * a) - (exp(a) / (2 * a) - sqrt(pi) * erfi(sqrt(a)) / (4 * a * sqrt(a)))
    A = findroot(lambda a: uA(a) - mpf(1) / 4, mpf(1.5))
    uB = lambda a: quad(lambda r: quad(lambda s: s * s * exp(a * s * s), [0, r]) / (r * r), [0, 1])
    B = findroot(lambda a: uB(a) - mpf(1) / 4, mpf(1.5))
    out = {"A Green's function closed form with erfi": A, "B iterated quadrature Int_0^1 r^-2 Int_0^r s^2 f": B}
elif k == 9:  # unit square, -div((1 + a x) grad u) = 1, u = 0 on x=0,1, Neumann on y=0,1: Int u = 1/20
    IA = lambda a: ((1 + a) * log(1 + a) - a) / (a * a * log(1 + a)) - 1 / (2 * a)
    A = findroot(lambda a: IA(a) - mpf(1) / 20, mpf(2))
    def IB(a):   # 1-D BVP -( (1+ax) u')' = 1 by BS shooting on u'(0), then quad of the interpolant-free integral
        def end(s):
            return bs_solve(lambda x, w: [w[1], (-1 - a * w[1]) / (1 + a * x), w[0]], 0, [0, s, 0], 1)
        s = findroot(lambda s: end(s)[0], mpf("0.4"))
        return end(s)[2]
    B = findroot(lambda a: IB(a) - mpf(1) / 20, A + mpf(10) ** -6)
    out = {"A closed form ((1+a)log(1+a) - a)/(a^2 log(1+a)) - 1/(2a)": A,
           "B Bulirsch-Stoer 1-D shooting with the integral as an extra ODE component": B}
elif k == 10:  # unit ball Bratu -Lap u = lam e^u, u = 0: u(0) = 1  -> lam = xi1^2/e (isothermal sphere)
    x0 = mpf(10) ** -12
    # universal isothermal-sphere equation w'' + 2 w'/x = e^-w, w(0) = w'(0) = 0, w = x^2/6 - x^4/120 + O(x^6)
    f = odefun(lambda x, w: [w[1], exp(-w[0]) - 2 * w[1] / x], x0, [x0 ** 2 / 6 - x0 ** 4 / 120, x0 / 3 - x0 ** 3 / 30])
    xi = findroot(lambda x: f(x)[0] - 1, mpf(2.8))
    A = xi ** 2 / exp(1)
    def u1(lam):
        c = lam * exp(1)
        return bs_solve(lambda r, w: [w[1], -2 * w[1] / r - lam * exp(w[0])], x0, [1 - c * x0 ** 2 / 6 + c * c * x0 ** 4 / 120, -c * x0 / 3 + c * c * x0 ** 3 / 30], 1)[0]
    B = findroot(u1, A + mpf(10) ** -8)
    out = {"A isothermal sphere (odefun) w(xi) = 1, lam = xi^2/e": A, "B Bulirsch-Stoer radial shooting in lam": B}
for kk, vv in out.items():
    print(kk); print("   ", nstr(vv, 70))
print("secs", round(time.time() - t0))
