"""Experiment 187 -- SIAM 100-Digit Challenge, App. D Problem 12 (Python column).

The eps -> 0 limiting ODE on the unit circle (derivation in constrained_particle.m):
    th'' = -cos th + alpha sin th cos th / (2 sqrt(1 + alpha cos^2 th)),
    th(0) = pi/2, th'(0) = -1;  find alpha with first hit of sin th = 0 at t = 10.

machine: scipy solve_ivp (DOP853, rtol = atol = 1e-12) with a terminal event
         on sin th = 0 (no hit by t = 100 -> 100), brentq on [10, 11] -- the
         same Brent-on-hitting-time statement as the WL row.
wp50:    mpmath has no event location, so the idiomatic statement is fixed-time
         shooting: odefun (Taylor) for th(10), findroot (secant) on
         sin th(10; alpha) = 0, started at 10.6234 / 10.62347 (both sides of
         the root, inside the window where the first hit is the only one).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

require(["scipy.integrate:solve_ivp", "scipy.optimize:brentq", "mpmath:odefun", "mpmath:findroot"])


def hit(a):
    def f(t, u):
        th, w = u
        c, s = math.cos(th), math.sin(th)
        return [w, -c + a * s * c / (2 * math.sqrt(1 + a * c * c))]
    ev = lambda t, u: math.sin(u[0])
    ev.terminal = True
    ev.direction = -1
    sol = solve_ivp(f, [0, 100], [math.pi / 2, -1.0], method="DOP853",
                    rtol=1e-12, atol=1e-12, events=ev)
    return sol.t_events[0][0] if len(sol.t_events[0]) else 100.0


vrun("d12.machine", lambda: brentq(lambda a: hit(a) - 10, 10, 11, xtol=1e-14), 10)


def y10(a):
    sol = mpmath.odefun(lambda t, u: [u[1], -mpmath.cos(u[0]) + a * mpmath.sin(u[0]) * mpmath.cos(u[0])
                                      / (2 * mpmath.sqrt(1 + a * mpmath.cos(u[0]) ** 2))],
                        0, [mpmath.pi / 2, mpmath.mpf(-1)])
    return mpmath.sin(sol(10)[0])


vrun("d12.wp50", lambda: mpmath.findroot(y10, (mpmath.mpf("10.6234"), mpmath.mpf("10.62347"))),
     48, cap=300, dps=60)
