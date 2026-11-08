"""Experiment 251 -- SIAM 100-digit SEALED HOLDOUT B, family D16 (Python column).

scipy has no PDE solver, so the idiomatic scipy statement is the method of
lines written by hand in the usual way (as exp 190): second-order finite
differences on a uniform grid of 600 intervals (the unit square: 100 x 100),
boundary conditions by ghost points, the initial data sampled at the nodes
(the mean of the one-sided values at a jump), solve_ivp (BDF, sparse
Jacobian, rtol = 1e-10) with a terminal EVENT where the discrete functional
(node value; trapezoid rule for Int u and Int u^2; the polygonal graph length
Sum sqrt(h^2 + du^2) for the arc length) reaches the target.  The same PDE as
the WL row, but the semi-discretisation is chosen here, not by a library, so
its truncation error (O(h^2)) is part of the scipy answer.  mpmath has no PDE
solver: every wp50 row is a SKIP.  See heat_functionals.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp

require(["scipy.integrate:solve_ivp", "scipy.sparse:diags"])

M = 600


def sample(f, x, eps=1e-12):
    return 0.5 * (f(x - eps) + f(x + eps))


def heat_time(u0, bc, func, target, tmax, L=1.0, src=0.0):
    h = L / M
    x = np.linspace(0.0, L, M + 1)
    if bc == "D":
        sl = slice(1, M); n = M - 1
        A = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n), format="lil")
    elif bc == "N":                                   # u_x = 0 at both ends (ghost points)
        sl = slice(0, M + 1); n = M + 1
        A = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n), format="lil")
        A[0, 1] = 2.0; A[n - 1, n - 2] = 2.0
    else:                                             # u(0) = 0, u_x(1) + u(1) = 0
        sl = slice(1, M + 1); n = M
        A = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n), format="lil")
        A[n - 1, n - 2] = 2.0; A[n - 1, n - 1] = -2.0 - 2.0 * h
    A = (A.tocsr()) / h ** 2

    def full(v):
        u = np.zeros(M + 1); u[sl] = v
        return u

    def ev(t, v):
        return func(full(v), x, h) - target
    ev.terminal = True
    v0 = sample(u0, x)[sl]
    sol = solve_ivp(lambda t, v: A @ v + src, [0, tmax], v0, method="BDF", jac=A,
                    rtol=1e-10, atol=1e-12, events=ev)
    return sol.t_events[0][0]


def point(x0):
    return lambda u, x, h: u[int(round(x0 / h))]


trap = lambda u, x, h: h * (u.sum() - 0.5 * (u[0] + u[-1]))
trap2 = lambda u, x, h: h * ((u ** 2).sum() - 0.5 * (u[0] ** 2 + u[-1] ** 2))
arc = lambda u, x, h: np.sqrt(h * h + np.diff(u) ** 2).sum()


def heat_time_2d(target, tmax, n=99):
    h = 1.0 / (n + 1)
    d2 = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n)) / h ** 2
    eye = sp.identity(n)
    lap = (sp.kron(d2, eye) + sp.kron(eye, d2)).tocsr()
    c = (n // 2) * n + n // 2
    ev = lambda t, v: v[c] - target
    ev.terminal = True
    sol = solve_ivp(lambda t, v: lap @ v, [0, tmax], np.ones(n * n), method="BDF", jac=lap,
                    rtol=1e-10, atol=1e-12, events=ev)
    return sol.t_events[0][0]


one = lambda x: np.ones_like(x)
haar_psi = lambda x: np.where((x >= 0) & (x < 0.5), 1.0, np.where((x >= 0.5) & (x < 1), -1.0, 0.0))
box3 = lambda x: np.where((x >= 1 / 3) & (x < 2 / 3), 1.0, 0.0)
CASES = [
    ("bd16.1", lambda: heat_time(one, "D", point(0.5), 0.5, 0.5)),
    ("bd16.2", lambda: heat_time(lambda x: 1 - np.abs(2 * x - 1), "D", trap, 0.25, 0.5)),
    ("bd16.3", lambda: heat_time(haar_psi, "D", arc, 3.0, 0.1)),
    ("bd16.4", lambda: heat_time(lambda x: x, "D", arc, 2.0, 0.1)),
    ("bd16.5", lambda: heat_time(lambda x: np.where(x > 1 / 3, 1.0, 0.0), "N", point(0.0), 1 / 3, 0.5)),
    ("bd16.6", lambda: heat_time(lambda x: x * (np.pi - x), "D", trap2, 1.0, 3.0, L=np.pi)),
    ("bd16.7", lambda: heat_time(lambda x: 0 * x, "D", point(0.5), 0.1, 0.5, src=1.0)),
    ("bd16.8", lambda: heat_time(one, "R", point(1.0), 0.25, 1.0)),
    ("bd16.9", lambda: heat_time_2d(0.5, 0.2)),
    ("bd16.10", lambda: heat_time(box3, "D", arc, 2.8, 0.02)),
]
for lab, f in CASES:
    vrun(lab + " machine", f, 10, cap=90)
    skip(lab + " wp50", "NDSolve (PDE)")
