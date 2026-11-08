"""Reference for holdout-B family d12 (experiment 247): find a parameter such that
an ODE EVENT (first hitting time) happens at a prescribed time.

Two independent methods per problem:
  A  closed form / quadrature formula for the event time + findroot, or (no
     closed form) mpmath.odefun (Taylor) fixed-time shooting + findroot;
  B  mpmath.odefun fixed-time shooting (when A is analytic), or the
     Gragg-Bulirsch-Stoer integrator of bhold_ode_util.py, fixed-time shooting.
"Fixed-time shooting": the event function evaluated at the prescribed time is
zero, e.g. th'(2; v) = 0; valid because each event function is monotone
through its first zero on the bracket (checked by a scipy scan when written).
Usage: bd12_ref.py k [dps]
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import (mp, mpf, sin, cos, sqrt, atan, acosh, atanh, ellipk, ellipf, asin, exp, log, quad,
                    odefun, findroot, nstr, pi, tan)
from bhold_ode_util import bs_solve

k = int(sys.argv[1]); mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 75


def at(rhs, y0, T, integ):
    return odefun(rhs, 0, y0)(T) if integ == "taylor" else bs_solve(rhs, 0, y0, T)


def fixed(rhs, y0f, T, g, x0, integ):
    return findroot(lambda p: g(p, at(rhs(p), y0f(p), T, integ)), mpf(x0))


out = {}
t0 = time.time()
if k == 1:   # pendulum th(0)=0, th'(0)=v; first th' = 0 at t = 2
    m = findroot(lambda m: ellipk(m) - 2, mpf("0.66")); A = 2 * sqrt(m)
    B = fixed(lambda v: (lambda t, u: [u[1], -sin(u[0])]), lambda v: [mpf(0), v], 2,
              lambda v, y: y[1], "1.62", "taylor")
    out = {"A v = 2 sqrt(m), K(m) = 2": A, "B odefun th'(2; v) = 0": B}
elif k == 2:  # vertical quadratic drag c = 1/10: landing time = 3
    c = mpf(1) / 10
    A = findroot(lambda v: atan(v * sqrt(c)) / sqrt(c) + acosh(sqrt(1 + c * v * v)) / sqrt(c) - 3, mpf("1.6"))
    B = findroot(lambda v: quad(lambda w: 1 / (1 + c * w * w), [0, v])
                 + quad(lambda w: 1 / (1 - c * w * w), [0, v / sqrt(1 + c * v * v)]) - 3, mpf("1.6"))
    out = {"A closed form atan + acosh": A, "B quad of the velocity-integral times (ascent + descent)": B}
elif k == 3:  # van der Pol x(0)=2, x'(0)=0: first x = 0 at t = 2
    rhs = lambda mu: (lambda t, u: [u[1], mu * (1 - u[0] ** 2) * u[1] - u[0]])
    A = fixed(rhs, lambda mu: [mpf(2), mpf(0)], 2, lambda mu, y: y[0], "0.8", "taylor")
    B = fixed(rhs, lambda mu: [mpf(2), mpf(0)], 2, lambda mu, y: y[0], "0.8", "bs")
    out = {"A odefun x(2; mu) = 0": A, "B Bulirsch-Stoer x(2; mu) = 0": B}
elif k == 4:  # Lotka-Volterra x(0)=a, y(0)=1: first x = 1 at t = 1
    rhs = lambda a: (lambda t, u: [u[0] * (1 - u[1]), u[1] * (u[0] - 1)])
    A = fixed(rhs, lambda a: [a, mpf(1)], 1, lambda a, y: y[0] - 1, "2.6", "taylor")
    B = fixed(rhs, lambda a: [a, mpf(1)], 1, lambda a, y: y[0] - 1, "2.6", "bs")
    out = {"A odefun x(1; a) = 1": A, "B Bulirsch-Stoer x(1; a) = 1": B}
elif k == 5:  # Kepler from (1,0) with velocity (0,v): first x = 0 at t = 13/10
    def tcross(v):
        a = 1 / (2 - v * v); e = 1 - v * v          # start = apoapsis (v < 1)
        nu = 3 * pi / 2
        E = 2 * pi - 2 * atan(sqrt((1 - e) / (1 + e)) / tan(pi / 4))   # tan(nu/2) = -1, E in (pi, 2pi)
        return sqrt(a ** 3) * ((E - e * sin(E)) - pi)
    A = findroot(lambda v: tcross(v) - mpf(13) / 10, mpf("0.7"))
    rhs = lambda v: (lambda t, u: [u[2], u[3], -u[0] / (u[0] ** 2 + u[1] ** 2) ** mpf(1.5),
                                   -u[1] / (u[0] ** 2 + u[1] ** 2) ** mpf(1.5)])
    B = fixed(rhs, lambda v: [mpf(1), mpf(0), mpf(0), v], mpf(13) / 10, lambda v, y: y[0], "0.7", "taylor")
    out = {"A Kepler equation from apoapsis (E - e sin E)": A, "B odefun x(13/10; v) = 0": B}
elif k == 6:  # damped oscillator x(0)=1, x'(0)=0: first x = 0 at t = 2
    def g(z):
        w = sqrt(1 - z * z); return 2 * w - (pi - atan(w / z))
    A = findroot(g, mpf("0.5"))
    B = fixed(lambda z: (lambda t, u: [u[1], -2 * z * u[1] - u[0]]), lambda z: [mpf(1), mpf(0)], 2,
              lambda z, y: y[0], "0.5", "taylor")
    out = {"A closed form: 2 w = pi - atan(w/zeta)": A, "B odefun x(2; zeta) = 0": B}
elif k == 7:  # Riccati y' = alpha + y^2, y(0) = 0: first y = 1 at t = 1
    A = findroot(lambda a: atan(1 / sqrt(a)) / sqrt(a) - 1, mpf("0.6"))
    B = fixed(lambda a: (lambda t, u: [a + u[0] ** 2]), lambda a: [mpf(0)], 1, lambda a, y: y[0] - 1,
              "0.6", "taylor")
    out = {"A closed form atan(1/sqrt a)/sqrt a = 1": A, "B odefun y(1; alpha) = 1": B}
elif k == 8:  # elastic pendulum, gravity 1, rest length 1, start (0.66,-0.88) at rest: first x=0 at t=2
    def rhs(kk):
        def f(t, u):
            x, y, vx, vy = u; r = sqrt(x * x + y * y); s = kk * (1 - 1 / r)
            return [vx, vy, -s * x, -s * y - 1]
        return f
    y0 = lambda kk: [mpf(66) / 100, mpf(-88) / 100, mpf(0), mpf(0)]
    A = fixed(rhs, y0, 2, lambda kk, y: y[0], "3.97", "taylor")
    B = fixed(rhs, y0, 2, lambda kk, y: y[0], "3.97", "bs")
    out = {"A odefun x(2; k) = 0": A, "B Bulirsch-Stoer x(2; k) = 0": B}
elif k == 9:  # SIR S(0)=0.99, I(0)=0.01: infection peak (beta S = 1) at t = 5
    S0, I0 = mpf(99) / 100, mpf(1) / 100
    def tpeak(b):
        I = lambda S: S0 + I0 - S + log(S / S0) / b
        return quad(lambda S: 1 / (b * S * I(S)), [1 / b, S0])
    A = findroot(lambda b: tpeak(b) - 5, mpf("2.0"))
    B = fixed(lambda b: (lambda t, u: [-b * u[0] * u[1], b * u[0] * u[1] - u[1]]), lambda b: [S0, I0], 5,
              lambda b, y: b * y[0] - 1, "2.0", "taylor")
    out = {"A quadrature t = Int dS/(beta S I(S)) on the SIR invariant": A, "B odefun beta S(5) = 1": B}
elif k == 10:  # pendulum th(0)=0, th'(0)=v: first th = 3 at t = 4 (just past the fold v ~ 2 sin(3/2))
    v0 = 2 * sin(mpf(3) / 2)
    A = findroot(lambda v: ellipf(asin(sin(mpf(3) / 2) / (v / 2)), (v / 2) ** 2) - 4,
                 (v0 + mpf(10) ** -12, v0 + mpf(10) ** -3), solver="anderson")
    B = fixed(lambda v: (lambda t, u: [u[1], -sin(u[0])]), lambda v: [mpf(0), v], 4,
              lambda v, y: y[0] - 3, nstr(A, 20), "taylor")
    out = {"A incomplete elliptic F(asin(sin(3/2)/k) | k^2) = 4, k = v/2": A, "B odefun th(4; v) = 3": B}
for kk, vv in out.items():
    print(kk); print("   ", nstr(vv, 70))
print("secs", round(time.time() - t0))
