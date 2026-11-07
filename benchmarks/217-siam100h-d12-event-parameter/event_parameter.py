"""Experiment 217 -- SIAM 100-Digit Challenge, holdout A, family d12 (Python column).

The ten problems of event_parameter.m: find a such that the first event of an
ODE fires at a given time T.
machine: scipy solve_ivp (DOP853, rtol = atol = 1e-12) with a terminal event
         (no event by t = 100 -> 100), brentq on the bracket of the .m file --
         the same Brent-on-hitting-time statement as the WL row.
wp50:    mpmath has no event location, so the idiomatic statement is
         fixed-time shooting (as in experiment 187): odefun (Taylor) for the
         state at T, findroot (Anderson, bracketing) on the event function at
         T over the same bracket (each bracket holds one sign change; checked).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

require(["scipy.integrate:solve_ivp", "scipy.optimize:brentq", "mpmath:odefun", "mpmath:findroot"])

# (rhs(t, u, a, lib), u0, event(u) as a function of the state, T, bracket)
# lib is math (machine) or mpmath (wp50); event is a level crossing in u.
def p9(t, u, a, m):
    c, s = m.cos(u[0]), m.sin(u[0])
    return [u[1], -c + a * s * c / (2 * m.sqrt(1 + a * c * c))]

P = {
    1: (lambda t, u, a, m: [u[1], -1 - a * u[1] ** 2], ["0", "2"], lambda u: u[1], "1", ("0.5", "3")),
    2: (lambda t, u, a, m: [u[1], -a * m.sin(u[0])], ["0", "1"], lambda u: u[1], "2", ("0.4", "2")),
    3: (lambda t, u, a, m: [u[1], -u[0] - a * u[0] ** 3], ["1", "0"], lambda u: u[0], "1.2", ("0", "3")),
    4: (lambda t, u, a, m: [u[1], -a * u[1] - u[0]], ["1", "0"], lambda u: u[0], "2", ("0.1", "1.5")),
    5: (lambda t, u, a, m: [a + u[0] ** 2], ["0"], lambda u: u[0] - 1, "0.5", ("1", "4")),
    6: (lambda t, u, a, m: [u[1], -a / u[0] ** 2], ["1", "0"], lambda u: u[0] - 0.5, "1", ("0.5", "1")),
    7: (lambda t, u, a, m: [u[1], -a * t * u[0]], ["1", "0"], lambda u: u[0], "1", ("3", "15")),
    8: (lambda t, u, a, m: [u[1], -(a + m.cos(t)) * u[0]], ["1", "0"], lambda u: u[0], "1", ("0.5", "3")),
    9: (p9, None, lambda u: math.sin(u[0]) if isinstance(u[0], float) else mpmath.sin(u[0]), "5", ("10.57", "10.58")),
    10: (lambda t, u, a, m: [u[1], a * (1 - u[0] ** 2) * u[1] - u[0]], ["2", "0"], lambda u: u[0], "2", ("0.2", "2")),
}


def u0_of(k, m):
    if k == 9:
        return [m.pi / 2, -1 if m is math else mpmath.mpf(-1)]
    return [float(v) if m is math else mpmath.mpf(v) for v in P[k][1]]


def hit(k, a):
    f, _, ev, _, _ = P[k]
    e = lambda t, u: ev(u)
    e.terminal = True
    sol = solve_ivp(lambda t, u: f(t, u, a, math), [0, 100], u0_of(k, math), method="DOP853",
                    rtol=1e-12, atol=1e-12, events=e)
    return sol.t_events[0][0] if len(sol.t_events[0]) else 100.0


def mp_case(k):
    f, _, ev, T, (lo, hi) = P[k]
    T = mpmath.mpf(T)
    def g(a):
        u = mpmath.odefun(lambda t, u: f(t, u, a, mpmath), 0, u0_of(k, mpmath))(T)
        if k == 6:
            return u[0] - mpmath.mpf(1) / 2
        return ev(u)
    return mpmath.findroot(g, (mpmath.mpf(lo), mpmath.mpf(hi)), solver="anderson")


for k in range(1, 11):
    T, (lo, hi) = float(P[k][3]), P[k][4]
    vrun("hd12.%d machine" % k, lambda: brentq(lambda a: hit(k, a) - T, float(lo), float(hi), xtol=1e-15), 10,
         cap=60)
    vrun("hd12.%d wp50" % k, lambda: mp_case(k), 48, cap=90, dps=60)
