"""Experiment 188 -- SIAM 100-Digit Challenge, App. D Problem 13 (Python column).

machine: scipy solve_ivp (DOP853, rtol = atol = 1e-12) with a terminal event
         on |u|^2 = 1 (never reached by t = 100 -> 100), brentq on z0 in
         [0.0064, 0.0065] for "first reaches |u| = 1 at t = 50" -- the same
         statement as the WL row (equivalent to |u(50)| = 1; see the .m header).
wp50:    mpmath has no event location, so the statement is the literal one:
         odefun (Taylor) for u(50), findroot (secant) on |u(50)| - 1, started at
         0.00648239617962 / 0.00648239617963 (|u(50)| is only defined on a
         relative window ~1e-12 wide: outside it the orbit has decayed or blown up).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

require(["scipy.integrate:solve_ivp", "scipy.optimize:brentq", "mpmath:odefun", "mpmath:findroot"])


def f(t, u):
    x, y, z = u
    n = math.sqrt(x * x + y * y + z * z)
    return [-x + 10 * y + n * (-0.7 * y - 0.03 * z),
            -y + 10 * z + n * (0.7 * y - 0.1 * z),
            -z + n * (0.03 * x + 0.1 * y)]


def reach(z0):
    ev = lambda t, u: u[0] * u[0] + u[1] * u[1] + u[2] * u[2] - 1
    ev.terminal = True
    sol = solve_ivp(f, [0, 100], [0.0, 0.0, z0], method="DOP853",
                    rtol=1e-12, atol=1e-12, events=ev)
    return sol.t_events[0][0] if len(sol.t_events[0]) else 100.0


vrun("d13.machine", lambda: brentq(lambda z0: reach(z0) - 50, 0.0064, 0.0065, xtol=1e-17), 10)


def norm50(z0):
    a7, a03, a1 = mpmath.mpf(7) / 10, mpmath.mpf(3) / 100, mpmath.mpf(1) / 10

    def F(t, u):
        x, y, z = u
        n = mpmath.sqrt(x * x + y * y + z * z)
        return [-x + 10 * y + n * (-a7 * y - a03 * z),
                -y + 10 * z + n * (a7 * y - a1 * z),
                -z + n * (a03 * x + a1 * y)]
    u = mpmath.odefun(F, 0, [mpmath.mpf(0), mpmath.mpf(0), z0])(50)
    return mpmath.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2) - 1


vrun("d13.wp50", lambda: mpmath.findroot(norm50, (mpmath.mpf("0.00648239617962"),
                                                  mpmath.mpf("0.00648239617963"))),
     48, cap=400, dps=60)
