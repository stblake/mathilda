"""Experiment 218 -- SIAM 100-Digit Challenge, holdout A, family d13 (Python column).

The ten problems of initial_value_shooting.m: find the initial value for
which the solution meets a condition at the final time T.
machine: scipy solve_ivp (DOP853, rtol = atol = 1e-12) for the state at T,
         brentq on the bracket of the .m file.  hd13.9 as in experiment 188:
         a terminal event on |u|^2 = 1 (never -> 100), brentq on reach - 30.
wp50:    mpmath odefun (Taylor) for the state at T, findroot (Anderson,
         bracketing) over the same bracket.  hd13.9: mpmath has no events, so
         the literal statement |u(30)| - 1 = 0 by secant from two starts in
         the window where |u(30)| is finite (as experiment 188 does).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

require(["scipy.integrate:solve_ivp", "scipy.optimize:brentq", "mpmath:odefun", "mpmath:findroot"])


def f9(t, u, m):
    x, y, z = u
    n = m.sqrt(x * x + y * y + z * z)
    c7, c03, c1 = (0.7, 0.03, 0.1) if m is math else (mpmath.mpf(7) / 10, mpmath.mpf(3) / 100, mpmath.mpf(1) / 10)
    return [-x + 10 * y + n * (-c7 * y - c03 * z), -y + 10 * z + n * (c7 * y - c1 * z), -z + n * (c03 * x + c1 * y)]


def kep(t, u, m):
    r3 = (u[0] ** 2 + u[1] ** 2) ** 1.5 if m is math else (u[0] ** 2 + u[1] ** 2) ** mpmath.mpf(1.5)
    return [u[2], u[3], -u[0] / r3, -u[1] / r3]

# k: (rhs(t, u, m), u0(p) as strings/p, T, residual(state, m), bracket)
P = {
    1: (lambda t, u, m: [u[0] * (1 - u[0])], lambda p: [p], 10, lambda u, m: u[0] - 0.5, (1e-5, 1e-4)),
    2: (lambda t, u, m: [-u[0] + 10 * u[1], -u[1] + 10 * u[2], -u[2]], lambda p: [0, 0, p], 5,
        lambda u, m: m.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) - 1, (0.05, 0.5)),
    3: (lambda t, u, m: [u[0] ** 2 - u[0]], lambda p: [p], 3, lambda u, m: u[0] - 2, (1, 1.04)),
    4: (lambda t, u, m: [u[1], -m.sin(u[0])], lambda p: [0, p], 2, lambda u, m: u[0] - 1, (0.5, 1.5)),
    5: (lambda t, u, m: [u[1], -m.exp(u[0])], lambda p: [0, p], 1, lambda u, m: u[0], (0.2, 1)),
    6: (lambda t, u, m: [u[1], u[2], -u[0] * u[2] / 2], lambda p: [0, 0, p], 10, lambda u, m: u[1] - 1, (0.2, 0.5)),
    7: (kep, lambda p: [1, 0, 0, p], 2, lambda u, m: m.sqrt(u[0] ** 2 + u[1] ** 2) - 1.5, (1.1, 1.3)),
    8: (lambda t, u, m: [t * t + u[0] ** 2], lambda p: [p], 1, lambda u, m: u[0] - 2, (0, 0.9)),
    9: (f9, lambda p: [0, 0, p], 30, lambda u, m: m.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) - 1, (0.0064, 0.0065)),
    10: (lambda t, u, m: [u[1], -u[0] - u[0] ** 3], lambda p: [p, 0], 10, lambda u, m: u[0], (0.3, 0.6)),
}


def machine(k):
    f, u0, T, res, (lo, hi) = P[k]
    if k == 9:
        def reach(p):
            ev = lambda t, u: u[0] ** 2 + u[1] ** 2 + u[2] ** 2 - 1
            ev.terminal = True
            s = solve_ivp(lambda t, u: f(t, u, math), [0, 100], [0.0, 0.0, p], method="DOP853",
                          rtol=1e-12, atol=1e-12, events=ev)
            return s.t_events[0][0] if len(s.t_events[0]) else 100.0
        return brentq(lambda p: reach(p) - 30, lo, hi, xtol=1e-17)
    def r(p):
        s = solve_ivp(lambda t, u: f(t, u, math), [0, T], [float(v) for v in u0(p)], method="DOP853",
                      rtol=1e-12, atol=1e-12)
        return res(s.y[:, -1], math)
    return brentq(r, lo, hi, xtol=1e-17)


def wp(k):
    f, u0, T, res, (lo, hi) = P[k]
    def r(p):
        u = mpmath.odefun(lambda t, u: f(t, u, mpmath), 0, [mpmath.mpf(v) for v in u0(p)])(T)
        if k == 7:
            return mpmath.sqrt(u[0] ** 2 + u[1] ** 2) - mpmath.mpf(3) / 2
        if k == 1:
            return u[0] - mpmath.mpf(1) / 2
        return res(u, mpmath)
    if k == 9:
        return mpmath.findroot(r, (mpmath.mpf("0.006483984"), mpmath.mpf("0.006483985")))
    return mpmath.findroot(r, (mpmath.mpf(lo), mpmath.mpf(hi)), solver="anderson")


for k in range(1, 11):
    vrun("hd13.%d machine" % k, lambda: machine(k), 10, cap=60)
    vrun("hd13.%d wp50" % k, lambda: wp(k), 48, cap=90, dps=60)
