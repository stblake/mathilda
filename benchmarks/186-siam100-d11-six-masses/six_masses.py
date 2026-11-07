"""Experiment 186 -- SIAM 100-Digit Challenge, Appendix D Problem 11 (Python column).

The shooting statement of six_masses.m: integrate r'' = Sum (m_i - r)/|m_i - r|^3
to t = T and solve {x(T), y(T)} = {4, 1} for (theta, T) from the same start
(-0.27, 1.92).  machine: scipy.integrate.solve_ivp (DOP853, rtol = atol =
1e-13) inside scipy.optimize.fsolve; wp50: mpmath.odefun (Taylor) inside
mpmath.findroot (multidimensional Newton) at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import numpy as np
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

require(["scipy.integrate:solve_ivp", "scipy.optimize:fsolve", "mpmath:odefun", "mpmath:findroot"])

M = np.array([[2, -1], [2, 0], [2, 1], [3, -1], [3, 0], [3, 1]], float)
def rhs(t, u):
    d = M - u[:2]
    a = (d / ((d * d).sum(1) ** 1.5)[:, None]).sum(0)
    return [u[2], u[3], a[0], a[1]]
def miss(z):
    th, T = z
    s = solve_ivp(rhs, [0, T], [0, 0, np.cos(th), np.sin(th)], method="DOP853", rtol=1e-13, atol=1e-13)
    return [s.y[0, -1] - 4, s.y[1, -1] - 1]
vrun("d11.machine", lambda: fsolve(miss, [-0.27, 1.92], xtol=1e-13)[1], 10)

def rhs_mp(t, u):
    x, y, vx, vy = u
    ax = ay = mpmath.mpf(0)
    for (a, b) in ((2, -1), (2, 0), (2, 1), (3, -1), (3, 0), (3, 1)):
        dx, dy = a - x, b - y
        r3 = (dx * dx + dy * dy) ** mpmath.mpf(1.5)
        ax += dx / r3; ay += dy / r3
    return [vx, vy, ax, ay]
def wp50():
    cache = {}
    def F(th, T):
        if th not in cache:
            cache[th] = mpmath.odefun(rhs_mp, 0, [0, 0, mpmath.cos(th), mpmath.sin(th)])
        u = cache[th](T)
        return [u[0] - 4, u[1] - 1]
    return mpmath.findroot(F, (mpmath.mpf("-0.27"), mpmath.mpf("1.92")))[1]
vrun("d11.wp50", wp50, 48, dps=60)
