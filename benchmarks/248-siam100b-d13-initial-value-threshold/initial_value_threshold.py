"""Experiment 248 -- SIAM 100-Digit Challenge holdout B for App. D Problem 13 (Python column).

The ten "find the initial value" problems of initial_value_threshold.m.
machine: event problems -- scipy solve_ivp (DOP853, rtol = atol = 1e-12) with a
         terminal event (none by t = 100 -> 100) and brentq on the machine
         bracket; fixed-time problems -- solve_ivp to T and scipy.optimize.newton
         (secant) from the same two starts as the WL row.
wp50:    mpmath has no event location, so every row is the fixed-time
         statement: odefun (Taylor) for the state at T, findroot (secant) on
         g(u(T)) = 0 from the narrow wp50 start pair.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from mpmath import mpf
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, newton

require(["scipy.integrate:solve_ivp", "scipy.optimize:brentq", "scipy.optimize:newton",
         "mpmath:odefun", "mpmath:findroot"])


def spec(M):
    sin, sq, pi = (mpmath.sin, mpmath.sqrt, mpmath.pi) if M else (math.sin, math.sqrt, math.pi)
    c = mpf if M else float
    a7, a03, a1 = (c(7) / 10, c(3) / 100, c(1) / 10)
    def model(t, u):
        x, y, z = u; n = sq(x * x + y * y + z * z)
        return [-x + 10 * y + n * (-a7 * y - a03 * z), -y + 10 * z + n * (a7 * y - a1 * z), -z + n * (a03 * x + a1 * y)]
    # k: (kind, rhs, y0(p), g(u) residual/event, T, machine bracket/starts, wp starts)
    return {
        1: ("ev", lambda t, u: [u[0] ** 2 - u[0]], lambda p: [p], lambda u: u[0] - c(1) / 2, 10,
            (0.99, 0.999999), ("0.99995", "0.99996")),
        2: ("fx", lambda t, u: [-u[0] + 10 * u[1], -u[1] + 10 * u[2], -u[2]], lambda p: [c(0), c(0), p],
            lambda u: sq(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) - 1, 10, (4, 5), ("4.4", "4.5")),
        3: ("fx", lambda t, u: [u[0] ** 2 - t], lambda p: [p], lambda u: u[0], 3, (0.727, 0.728),
            ("0.7276", "0.7278")),
        4: ("ev", model, lambda p: [c(0), c(0), p], lambda u: u[0] ** 2 + u[1] ** 2 + u[2] ** 2 - c(1) / 4, 40,
            (0.0064, 0.0065), ("0.006482442672312", "0.006482442672313")),
        5: ("ev", lambda t, u: [u[1], -sin(u[0]) - u[1] / 10], lambda p: [c(0), p], lambda u: u[0] - 4, 10,
            (2.2, 2.3), ("2.204", "2.2041")),
        6: ("ev", lambda t, u: [u[1], -1 / u[0] ** 2], lambda p: [c(1), p], lambda u: u[0] - 3, 5,
            (1.1, 1.3), ("1.15503", "1.15505")),
        7: ("ev", lambda t, u: [u[0] * (1 - u[0]) * (u[0] - c(3) / 10)], lambda p: [p], lambda u: u[0] - c(9) / 10,
            10, (0.31, 0.5), ("0.368", "0.369")),
        8: ("fx", lambda t, u: [u[1], u[2], -u[0] * u[2] / 2], lambda p: [c(0), c(0), p], lambda u: u[1] - 1, 10,
            (0.3, 0.4), ("0.332", "0.333")),
        9: ("ev", lambda t, u: [u[0] * (1 - u[0]) - (1 + sin(2 * pi * t)) / 10], lambda p: [p],
            lambda u: u[0] - c(1) / 2, 10, (0.12, 0.2), ("0.1288", "0.1289")),
        10: ("fx", lambda t, u: [u[1], u[2], -u[0] * u[2] - 1 + u[1] ** 2], lambda p: [c(0), c(0), p],
             lambda u: u[1] - 1, 10, (1.2, 1.25), ("1.2325", "1.2326")),
    }


def m_case(k):
    kind, rhs, y0, g, T, br, _ = spec(False)[k]
    if kind == "ev":
        def hit(p):
            ev = lambda t, u: g(u)
            ev.terminal = True
            s = solve_ivp(rhs, [0, 100], y0(p), method="DOP853", rtol=1e-12, atol=1e-12, events=ev)
            return s.t_events[0][0] if len(s.t_events[0]) else 100.0
        return brentq(lambda p: hit(p) - T, br[0], br[1], xtol=1e-17)
    def res(p):
        s = solve_ivp(rhs, [0, T], y0(p), method="DOP853", rtol=1e-12, atol=1e-12)
        return g(s.y[:, -1])
    return newton(res, br[0], x1=br[1], tol=1e-15)


def w_case(k):
    kind, rhs, y0, g, T, _, wb = spec(True)[k]
    F = lambda p: g(mpmath.odefun(rhs, 0, y0(p))(T))
    return mpmath.findroot(F, (mpf(wb[0]), mpf(wb[1])))


for k in range(1, 11):
    vrun("bd13.%d machine" % k, lambda k=k: m_case(k), 10)
for k in range(1, 11):
    vrun("bd13.%d wp50" % k, lambda k=k: w_case(k), 48, dps=60)
