"""Experiment 247 -- SIAM 100-Digit Challenge holdout B for App. D Problem 12 (Python column).

The ten "find p such that the first event happens at t = T*" problems of
event_time_root.m (see its header for the statements).
machine: scipy solve_ivp (DOP853, rtol = atol = 1e-12) with a terminal event
         (no event by t = 100 -> 100), brentq on the machine bracket -- the
         same Brent-on-event-time statement as the WL row.
wp50:    mpmath has no event location, so the idiomatic statement is fixed-time
         shooting: odefun (Taylor) for the state at T*, findroot (secant) on the
         event function there (e.g. th'(2; v) = 0), started at the two ends of
         the narrow wp50 bracket (as in exp. 187).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from mpmath import mpf
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

require(["scipy.integrate:solve_ivp", "scipy.optimize:brentq", "mpmath:odefun", "mpmath:findroot"])


def P(M):
    """Per problem: (rhs(p), y0(p), event(y), direction, target T, machine bracket, wp bracket)."""
    sin = mpmath.sin if M else math.sin
    sq = mpmath.sqrt if M else math.sqrt
    c = mpf if M else float
    return {
        1: (lambda p: lambda t, u: [u[1], -sin(u[0])], lambda p: [c(0), p], lambda u: u[1], -1, 2,
            (1, 1.9), ("1.60", "1.61")),
        2: (lambda p: lambda t, u: [u[1], -1 - u[1] * abs(u[1]) / 10], lambda p: [c(0), p], lambda u: u[0], -1, 3,
            (0.5, 5), ("1.58", "1.59")),
        3: (lambda p: lambda t, u: [u[1], p * (1 - u[0] ** 2) * u[1] - u[0]], lambda p: [c(2), c(0)],
            lambda u: u[0], -1, 2, (0, 2), ("0.75", "0.76")),
        4: (lambda p: lambda t, u: [u[0] * (1 - u[1]), u[1] * (u[0] - 1)], lambda p: [p, c(1)],
            lambda u: u[0] - 1, -1, 1, (1.1, 6), ("2.63", "2.64")),
        5: (lambda p: lambda t, u: [u[2], u[3], -u[0] / (u[0] ** 2 + u[1] ** 2) ** c(1.5),
                                    -u[1] / (u[0] ** 2 + u[1] ** 2) ** c(1.5)],
            lambda p: [c(1), c(0), c(0), p], lambda u: u[0], -1, c(13) / 10, (0.2, 0.95), ("0.57", "0.58")),
        6: (lambda p: lambda t, u: [u[1], -2 * p * u[1] - u[0]], lambda p: [c(1), c(0)], lambda u: u[0], -1, 2,
            (0, 0.9), ("0.31", "0.32")),
        7: (lambda p: lambda t, u: [p + u[0] ** 2], lambda p: [c(0)], lambda u: u[0] - 1, 1, 1,
            (0.1, 2), ("0.74", "0.75")),
        8: (lambda p: lambda t, u: [u[2], u[3], -p * (1 - 1 / sq(u[0] ** 2 + u[1] ** 2)) * u[0],
                                    -p * (1 - 1 / sq(u[0] ** 2 + u[1] ** 2)) * u[1] - 1],
            lambda p: [c(66) / 100, c(-88) / 100, c(0), c(0)], lambda u: u[0], -1, 2, (2, 5), ("3.98", "3.99")),
        9: (lambda p: lambda t, u: [-p * u[0] * u[1], p * u[0] * u[1] - u[1]], lambda p: [c(99) / 100, c(1) / 100],
            lambda p_u: None, -1, 5, (1.5, 3), ("1.81", "1.82")),
        10: (lambda p: lambda t, u: [u[1], -sin(u[0])], lambda p: [c(0), p], lambda u: u[0] - 3, 1, 4,
             (1.9, 2.1), ("1.99499", "1.995")),
    }


def event_fn(k, p, spec):
    if k == 9:
        return lambda t, u: p * u[0] - 1
    g = spec[2]
    return lambda t, u: g(u)


def m_case(k):
    spec = P(False)[k]
    rhs, y0, _, direc, T, br, _ = spec
    def hit(p):
        ev = event_fn(k, p, spec)
        ev.terminal = True
        ev.direction = direc
        s = solve_ivp(rhs(p), [0, 100], y0(p), method="DOP853", rtol=1e-12, atol=1e-12, events=ev)
        return s.t_events[0][0] if len(s.t_events[0]) else 100.0
    return brentq(lambda p: hit(p) - T, br[0], br[1], xtol=1e-14)


def w_case(k):
    spec = P(True)[k]
    rhs, y0, g, _, T, _, wb = spec
    def F(p):
        u = mpmath.odefun(rhs(p), 0, y0(p))(T)
        return p * u[0] - 1 if k == 9 else g(u)
    return mpmath.findroot(F, (mpf(wb[0]), mpf(wb[1])))


for k in range(1, 11):
    vrun("bd12.%d machine" % k, lambda k=k: m_case(k), 10)
for k in range(1, 11):
    vrun("bd12.%d wp50" % k, lambda k=k: w_case(k), 48, dps=60)
