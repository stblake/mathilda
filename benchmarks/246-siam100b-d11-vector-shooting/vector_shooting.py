"""Experiment 246 -- SIAM 100-Digit Challenge holdout B for App. D Problem 11 (Python column).

The shooting statements of vector_shooting.m (see its header): an ODE solve to
an unknown final time inside a vector root finder, from the same guesses.
machine: scipy.integrate.solve_ivp (DOP853, rtol = atol = 1e-13) inside
         scipy.optimize.fsolve (xtol 1e-13);
wp50:    mpmath.odefun (Taylor) inside mpmath.findroot (multidimensional
         Newton) at 60 digits.  odefun integrates forward only, so a Newton
         step to a negative time is a failure of that row, not of the script.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import numpy as np
import mpmath
from mpmath import mpf
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

require(["scipy.integrate:solve_ivp", "scipy.optimize:fsolve", "mpmath:odefun", "mpmath:findroot"])

MU = 12277471 / 10 ** 9


def acc_np(k, x, y, vx, vy):
    if k == 1: return -vx / 2, -1 - vy / 2
    if k == 2: return -x, -y
    if k == 3:
        r3 = (x * x + y * y) ** 1.5; return -x / r3, -y / r3
    if k == 7: return -x - 2 * x * y, -y - x * x + y * y
    if k == 8:
        ax = ay = 0.0
        for p in (-1.0, 1.0):
            dx, dy = p - x, -y; r3 = (dx * dx + dy * dy) ** 1.5; ax += dx / r3; ay += dy / r3
        return ax, ay


def acc_mp(k, x, y, vx, vy):
    if k == 1: return -vx / 2, -1 - vy / 2
    if k == 2: return -x, -y
    if k == 3:
        r2 = x * x + y * y; r3 = r2 * mpmath.sqrt(r2); return -x / r3, -y / r3
    if k == 7: return -x - 2 * x * y, -y - x * x + y * y
    if k == 8:
        ax = ay = mpf(0)
        for p in (-1, 1):
            dx, dy = p - x, -y; r2 = dx * dx + dy * dy; r3 = r2 * mpmath.sqrt(r2); ax += dx / r3; ay += dy / r3
        return ax, ay


PLANAR = {1: ((0, 0), 2, (1, 0.5), (0.45, 0.6)), 2: ((1, 0), 1, (0, 0.5), (2.47, 0.91)),
          3: ((1, 0), 1.1, (-0.5, 1), (1.71, 1.84)), 7: ((0, 0), 0.5, (0.3, 0.2), (0.5836, 0.8305)),
          8: ((0, -2), 0.8, (0.5, 1.5), (1.50, 2.75))}
def planar_mp():
    # built at call time, so the exact data are rounded at the case's 60 digits
    return {1: ((0, 0), 2, (1, mpf(1) / 2)), 2: ((1, 0), 1, (0, mpf(1) / 2)),
             3: ((1, 0), mpf(11) / 10, (mpf(-1) / 2, 1)), 7: ((0, 0), mpf(1) / 2, (mpf(3) / 10, mpf(1) / 5)),
             8: ((0, -2), mpf(4) / 5, (mpf(1) / 2, mpf(3) / 2))}


def ivp(f, T, y0):
    s = solve_ivp(f, [0, T], y0, method="DOP853", rtol=1e-13, atol=1e-13)
    return s.y[:, -1]


def odef(f, T, y0):
    return mpmath.odefun(f, 0, y0)(T)


def m_planar(k):
    (x0, y0), v, tgt, gs = PLANAR[k]
    rhs = lambda t, u: [u[2], u[3], *acc_np(k, *u)]
    def miss(z):
        th, T = z
        e = ivp(rhs, T, [x0, y0, v * np.cos(th), v * np.sin(th)])
        return [e[0] - tgt[0], e[1] - tgt[1]]
    return fsolve(miss, gs, xtol=1e-13)[1]


def w_planar(k):
    (x0, y0), v, tgt = planar_mp()[k]; gs = PLANAR[k][3]
    rhs = lambda t, u: [u[2], u[3], *acc_mp(k, *u)]
    def F(th, T):
        e = odef(rhs, T, [mpf(x0), mpf(y0), v * mpmath.cos(th), v * mpmath.sin(th)])
        return [e[0] - tgt[0], e[1] - tgt[1]]
    return mpmath.findroot(F, [mpf(gs[0]), mpf(gs[1])])[1]


# the non-planar problems: (rhs, initial state from the unknowns, final time index,
# residual, guesses, which unknown is reported)
def spec(k, M):
    sin, cos, exp = (np.sin, np.cos, np.exp) if not M else (mpmath.sin, mpmath.cos, mpmath.exp)
    c = (lambda a: a) if not M else mpf
    if k == 4:
        return (lambda t, u: [u[1], -u[2] * exp(u[0]), 0 * u[0]],
                lambda z: [c(0), z[1], z[0]], None, (2.7, 1.9), 0)
    if k == 5:
        return (lambda t, u: [u[1], u[2], u[3], 1 - 50 * (1 + t) * u[0]],
                lambda z: [c(0), c(0), z[0], z[1]], 1, (0.08, -0.5), 0)
    if k == 6:
        return (lambda t, u: [u[1], -sin(u[0])], lambda z: [c(0), z[0]], "T", (1.7, 2.1), 1)
    if k == 9:
        return (lambda t, u: [u[3], u[4], u[5], -u[0], -4 * u[1], -9 * u[2]],
                lambda z: [c(0)] * 3 + [cos(z[1]) * cos(z[0]), cos(z[1]) * sin(z[0]), sin(z[1])],
                "T", (0.6, 0.25, 0.38), 2)
    if k == 10:
        mu = MU if not M else mpf(12277471) / 10 ** 9
        m = 1 - mu
        sq = np.sqrt if not M else mpmath.sqrt
        def rhs(t, u):
            x, y, vx, vy = u
            d1 = (x + mu) ** 2 + y * y; d1 = d1 * sq(d1); d2 = (x - m) ** 2 + y * y; d2 = d2 * sq(d2)
            return [vx, vy, x + 2 * vy - m * (x + mu) / d1 - mu * (x - m) / d2,
                    y - 2 * vx - m * y / d1 - mu * y / d2]
        return (rhs, lambda z: [c(994) / 1000, c(0), c(0), z[0]], "T",
                (-2.0015851063790825, 8.532608280078981), 1)


def resid(k, z, e, M):
    half = 0.5 if not M else mpf(1) / 2
    if k == 4: return None
    if k == 5: return [e[0], e[1]]
    if k == 6: return [e[0] - 2, e[1]]
    if k == 9: return [e[0] - (0.3 if not M else mpf(3) / 10), e[1] - (0.2 if not M else mpf(1) / 5),
                       e[2] - (0.1 if not M else mpf(1) / 10)]
    if k == 10: return [e[1], e[2]]


def m_other(k):
    rhs, init, Tidx, gs, rep = spec(k, False)
    def miss(z):
        if k == 4:
            s = solve_ivp(rhs, [0, 1], init(z), method="DOP853", rtol=1e-13, atol=1e-13, t_eval=[0.5, 1])
            return [s.y[0, 1], s.y[0, 0] - 0.5]
        T = z[-1] if Tidx == "T" else Tidx
        return resid(k, z, ivp(rhs, T, init(z)), False)
    return fsolve(miss, gs, xtol=1e-13)[rep]


def w_other(k):
    rhs, init, Tidx, gs, rep = spec(k, True)
    def F(*z):
        if k == 4:
            f = mpmath.odefun(rhs, 0, init(z))
            return [f(1)[0], f(mpf(1) / 2)[0] - mpf(1) / 2]
        T = z[-1] if Tidx == "T" else Tidx
        return resid(k, z, odef(rhs, T, init(z)), True)
    return mpmath.findroot(F, [mpf(g) for g in gs])[rep]


for tag, d, dps in (("machine", 10, None), ("wp50", 48, 60)):
    for k in range(1, 11):
        if k in PLANAR:
            fn = (lambda k=k: m_planar(k)) if dps is None else (lambda k=k: w_planar(k))
        else:
            fn = (lambda k=k: m_other(k)) if dps is None else (lambda k=k: w_other(k))
        vrun("bd11.%d %s" % (k, tag), fn, d, dps=dps)
