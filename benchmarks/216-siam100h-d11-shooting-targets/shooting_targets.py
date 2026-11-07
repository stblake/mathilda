"""Experiment 216 -- SIAM 100-Digit Challenge, holdout A, family d11 (Python column).

The two-unknown shooting statement of shooting_targets.m: integrate the force
law from the start with speed v in direction theta to t = T and solve
{x(T), y(T)} = target for (theta, T) from the same starting guess.
machine: scipy.integrate.solve_ivp (DOP853, rtol = atol = 1e-13) inside
         scipy.optimize.fsolve (xtol 1e-13);
wp50:    mpmath.odefun (Taylor) inside mpmath.findroot (multidimensional
         Newton) at 60 digits, one odefun object cached per theta (as exp 186).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

require(["scipy.integrate:solve_ivp", "scipy.optimize:fsolve", "mpmath:odefun", "mpmath:findroot"])


def grav(ms):
    def acc(x, y, vx, vy, m):
        ax = ay = 0
        for (a, b) in ms:
            dx, dy = a - x, b - y
            r3 = (dx * dx + dy * dy) ** (1.5 if m is math else mpmath.mpf(1.5))
            ax += dx / r3
            ay += dy / r3
        return ax, ay
    return acc


def kepler(x, y, vx, vy, m):
    r3 = (x * x + y * y) ** (1.5 if m is math else mpmath.mpf(1.5))
    return -x / r3, -y / r3

# k: (acceleration(x, y, vx, vy, lib), start, speed, target, guess) -- speed/target as strings
C = {
    1: (lambda x, y, vx, vy, m: (0, -1), (0, 0), "2", ("1", "0.5"), ("0.6", "0.6")),
    2: (lambda x, y, vx, vy, m: (-x, -y), (1, 0), "1", ("0", "0.5"), ("2.5", "0.9")),
    3: (lambda x, y, vx, vy, m: (-x, -4 * y), (0, 0), "1", ("0.5", "0.25"), ("0.55", "0.63")),
    4: (lambda x, y, vx, vy, m: (-vx / 2, -vy / 2 - 1), (0, 0), "2", ("1", "0"), ("0.15", "0.58")),
    5: (kepler, (1, 0), "1", ("-0.5", "0.5"), ("1.89", "1.50")),
    6: (kepler, (1, 0), "1.2", ("0", "1.5"), ("1.54", "1.80")),
    7: (grav([(2, -1), (2, 1)]), (0, 0), "1", ("4", "0"), ("0.05", "2.9")),
    8: (grav([(1, 1), (2, -1), (3, 1)]), (0, 0), "1", ("4", "0"), ("-0.35", "2.88")),
    9: (lambda x, y, vx, vy, m: (vy, -vx), (0, 0), "1", ("1", "0.5"), ("1.06", "1.19")),
    10: (grav([(2, -1), (2, 0), (2, 1), (3, -1), (3, 0), (3, 1)]), (0, 0), "1", ("1", "2"), ("0.94", "2.03")),
}


def machine(k):
    acc, s0, v, tg, g = C[k]
    v = float(v)
    def rhs(t, u):
        a = acc(u[0], u[1], u[2], u[3], math)
        return [u[2], u[3], a[0], a[1]]
    def miss(z):
        th, T = z
        s = solve_ivp(rhs, [0, T], [s0[0], s0[1], v * math.cos(th), v * math.sin(th)], method="DOP853",
                      rtol=1e-13, atol=1e-13)
        return [s.y[0, -1] - float(tg[0]), s.y[1, -1] - float(tg[1])]
    return fsolve(miss, [float(g[0]), float(g[1])], xtol=1e-13)[1]


def wp(k):
    acc, s0, v, tg, g = C[k]
    v = mpmath.mpf(v)
    tg = [mpmath.mpf(c) for c in tg]
    def rhs(t, u):
        a = acc(u[0], u[1], u[2], u[3], mpmath)
        return [u[2], u[3], a[0], a[1]]
    cache = {}
    def F(th, T):
        if th not in cache:
            cache[th] = mpmath.odefun(rhs, 0, [mpmath.mpf(s0[0]), mpmath.mpf(s0[1]), v * mpmath.cos(th), v * mpmath.sin(th)])
        u = cache[th](T)
        return [u[0] - tg[0], u[1] - tg[1]]
    return mpmath.findroot(F, (mpmath.mpf(g[0]), mpmath.mpf(g[1])))[1]


for k in range(1, 11):
    vrun("hd11.%d machine" % k, lambda: machine(k), 10, cap=60)
    vrun("hd11.%d wp50" % k, lambda: wp(k), 48, cap=90, dps=60)
