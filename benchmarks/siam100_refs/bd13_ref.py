"""Reference for holdout-B family d13 (experiment 248): find an INITIAL value such
that the solution meets a condition at a prescribed later time (often through
the stable manifold of a saddle, so the map is very sensitive).

Two independent methods per problem:
  A  closed form (Airy functions, matrix exponential, Kepler's radial
     equation, partial fractions) + findroot, or else mpmath.odefun (Taylor)
     fixed-time shooting + findroot;
  B  mpmath.odefun fixed-time shooting (when A is analytic), or the
     Gragg-Bulirsch-Stoer integrator of bhold_ode_util.py.
Event problems ("first time g = 0 is t = T") are solved here in the
equivalent fixed-time form g(u(T)) = 0 (g is monotone through its first zero
on the bracket; see the experiment header).
Usage: bd13_ref.py k [dps]
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import (mp, mpf, sin, cos, sqrt, exp, log, acos, airyai, airybi, pi, odefun, findroot, nstr)
from bhold_ode_util import bs_solve

k = int(sys.argv[1]); mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 75


def at(rhs, y0, T, integ):
    return odefun(rhs, 0, y0)(T) if integ == "taylor" else bs_solve(rhs, 0, y0, T)


def fixed(rhs, y0f, T, g, x0, integ):
    """x0: a start (secant) or a (lo, hi) bracket (Anderson-Bjorck, stays inside)."""
    if isinstance(x0, tuple):
        return findroot(lambda p: g(p, at(rhs(p), y0f(p), T, integ)), (mpf(x0[0]), mpf(x0[1])), solver="anderson")
    return findroot(lambda p: g(p, at(rhs(p), y0f(p), T, integ)), mpf(x0))


out = {}
t0 = time.time()
if k == 1:   # y' = y^2 - y; first y = 1/2 at t = 10
    A = 1 / (1 + exp(-10))
    B = fixed(lambda p: (lambda t, u: [u[0] ** 2 - u[0]]), lambda p: [p], 10, lambda p, y: y[0] - mpf(1) / 2,
              ("0.9999", "0.99999"), "taylor")
    out = {"A 1/(1 + e^-10)": A, "B odefun y(10) = 1/2": B}
elif k == 2:  # u' = A u, A = -I + 10 (superdiagonal), u(0) = (0,0,z0): |u(10)| = 1
    A = exp(10) / sqrt(mpf(25010001))
    rhs = lambda p: (lambda t, u: [-u[0] + 10 * u[1], -u[1] + 10 * u[2], -u[2]])
    B = fixed(rhs, lambda p: [mpf(0), mpf(0), p], 10, lambda p, y: sqrt(y[0] ** 2 + y[1] ** 2 + y[2] ** 2) - 1,
              "4.4", "taylor")
    out = {"A e^10/sqrt(25010001) (exp of a Jordan block)": A, "B odefun |u(10)| = 1": B}
elif k == 3:  # Riccati y' = y^2 - t: y(3) = 0  -> y(0)
    r = -airyai(3, 1) / airybi(3, 1)
    A = -(airyai(0, 1) + r * airybi(0, 1)) / (airyai(0) + r * airybi(0))
    B = fixed(lambda p: (lambda t, u: [u[0] ** 2 - t]), lambda p: [p], 3, lambda p, y: y[0], ("0.7270", "0.7280"), "taylor")
    out = {"A Airy: y = -w'/w, w = Ai + r Bi, w'(3) = 0": A, "B odefun y(3) = 0": B}
elif k == 4:  # exp. 188's nonnormal model, |u| first = 1/2 at t = 40
    a7, a03, a1 = mpf(7) / 10, mpf(3) / 100, mpf(1) / 10
    def rhs(p):
        def F(t, u):
            x, y, z = u; n = sqrt(x * x + y * y + z * z)
            return [-x + 10 * y + n * (-a7 * y - a03 * z), -y + 10 * z + n * (a7 * y - a1 * z), -z + n * (a03 * x + a1 * y)]
        return F
    g = lambda p, y: y[0] ** 2 + y[1] ** 2 + y[2] ** 2 - mpf(1) / 4
    integ = sys.argv[3] if len(sys.argv) > 3 else "taylor"
    A = findroot(lambda p: g(p, at(rhs(p), [mpf(0), mpf(0), p], 40, integ)),
                 (mpf("0.0064824426723120"), mpf("0.0064824426723125")), solver="secant")
    out = {integ + " |u(40)|^2 = 1/4": A}
elif k == 5:  # damped pendulum th'' = -sin th - th'/10, th(0) = 0: first th = 4 at t = 10
    rhs = lambda v: (lambda t, u: [u[1], -sin(u[0]) - u[1] / 10])
    A = fixed(rhs, lambda v: [mpf(0), v], 10, lambda v, y: y[0] - 4, "2.2040396461", "taylor")
    B = fixed(rhs, lambda v: [mpf(0), v], 10, lambda v, y: y[0] - 4, "2.2040396461", "bs")
    out = {"A odefun th(10) = 4": A, "B Bulirsch-Stoer th(10) = 4": B}
elif k == 6:  # radial Kepler r'' = -1/r^2, r(0) = 1: first r = 3 at t = 5
    def tof(v):
        a = 1 / (2 - v * v); e0 = acos(1 - 1 / a); e1 = acos(1 - 3 / a)
        return a ** mpf(1.5) * ((e1 - sin(e1)) - (e0 - sin(e0)))
    A = findroot(lambda v: tof(v) - 5, mpf("1.15504"))
    B = fixed(lambda v: (lambda t, u: [u[1], -1 / u[0] ** 2]), lambda v: [mpf(1), v], 5, lambda v, y: y[0] - 3,
              "1.15504", "taylor")
    out = {"A radial Kepler r = a(1 - cos eta), t = a^(3/2)(eta - sin eta)": A, "B odefun r(5) = 3": B}
elif k == 7:  # Allee x' = x(1-x)(x - 3/10): first x = 9/10 at t = 10
    a = mpf(3) / 10
    F = lambda x: -log(x) / a - log(1 - x) / (1 - a) + log(x - a) / (a * (1 - a))
    A = findroot(lambda x0: F(mpf(9) / 10) - F(x0) - 10, mpf("0.3683"))
    B = fixed(lambda p: (lambda t, u: [u[0] * (1 - u[0]) * (u[0] - a)]), lambda p: [p], 10,
              lambda p, y: y[0] - mpf(9) / 10, "0.3683", "taylor")
    out = {"A partial fractions t = [ -ln x/a - ln(1-x)/(1-a) + ln(x-a)/(a(1-a)) ]": A, "B odefun x(10) = 9/10": B}
elif k == 8:  # Blasius f''' = -f f''/2 on [0,10]: f'(10) = 1 -> f''(0)
    rhs = lambda s: (lambda t, u: [u[1], u[2], -u[0] * u[2] / 2])
    A = fixed(rhs, lambda s: [mpf(0), mpf(0), s], 10, lambda s, y: y[1] - 1, "0.3320573372", "taylor")
    B = fixed(rhs, lambda s: [mpf(0), mpf(0), s], 10, lambda s, y: y[1] - 1, "0.3320573372", "bs")
    out = {"A odefun f'(10) = 1": A, "B Bulirsch-Stoer f'(10) = 1": B}
elif k == 9:  # x' = x(1-x) - (1/10)(1 + sin 2 pi t): first x = 1/2 at t = 10
    rhs = lambda p: (lambda t, u: [u[0] * (1 - u[0]) - (1 + sin(2 * pi * t)) / 10])
    A = fixed(rhs, lambda p: [p], 10, lambda p, y: y[0] - mpf(1) / 2, "0.12885281", "taylor")
    B = fixed(rhs, lambda p: [p], 10, lambda p, y: y[0] - mpf(1) / 2, "0.12885281", "bs")
    out = {"A odefun x(10) = 1/2": A, "B Bulirsch-Stoer x(10) = 1/2": B}
elif k == 10:  # Falkner-Skan beta = 1 (Hiemenz): f''' = -f f'' - 1 + f'^2 on [0,10], f'(10) = 1
    rhs = lambda s: (lambda t, u: [u[1], u[2], -u[0] * u[2] - 1 + u[1] ** 2])
    A = fixed(rhs, lambda s: [mpf(0), mpf(0), s], 10, lambda s, y: y[1] - 1, "1.2325876", "taylor")
    B = fixed(rhs, lambda s: [mpf(0), mpf(0), s], 10, lambda s, y: y[1] - 1, "1.2325876", "bs")
    out = {"A odefun f'(10) = 1": A, "B Bulirsch-Stoer f'(10) = 1": B}
for kk, vv in out.items():
    print(kk); print("   ", nstr(vv, 70))
print("secs", round(time.time() - t0))
