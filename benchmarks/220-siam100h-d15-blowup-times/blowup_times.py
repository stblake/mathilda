"""Experiment 220 -- SIAM 100-Digit Challenge, holdout A, family d15 (Python column).

machine: the ODE rows: scipy solve_ivp (DOP853, rtol = atol = 1e-12) on
         [0, tmax]; at the blow-up the solver gives up ("Required step size is
         less than spacing between numbers"), and that time -- like the end of
         NDSolve's domain -- is the answer.  The PDE rows (hd15.7, 8, 9): scipy
         has no PDE solver, so (as in experiment 190) the method of lines is
         written in the usual way: second-order finite differences
         (scipy.sparse; Neumann walls by the ghost-point rule), solve_ivp BDF
         with the sparse Jacobian, rtol = atol = 1e-10; hd15.9 then brentq on
         the dense output at the centre node.  The semi-discretisation is
         chosen here, not by a library, so its truncation error is part of the
         answer (exact for hd15.7/8, whose solution is spatially constant).
wp50:    mpmath has neither a blow-up-detecting ODE solver nor a PDE solver:
         SKIP.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import math
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

require(["scipy.integrate:solve_ivp", "scipy.sparse:kron", "scipy.optimize:brentq"])


def tinf(f, u0, tmax):
    s = solve_ivp(f, [0, tmax], u0, method="DOP853", rtol=1e-12, atol=1e-12)
    return s.t[-1] if s.status == -1 else "no blow-up"


ODE = {
    1: (lambda t, u: [1 + u[0] ** 2], [0.0], 3),
    2: (lambda t, u: [t * t + u[0] ** 2], [0.0], 3),
    3: (lambda t, u: [u[1], u[0] ** 2], [1.0, 0.0], 4),
    4: (lambda t, u: [math.exp(u[0]) - u[0]], [0.0], 3),
    5: (lambda t, u: [u[0] ** 2 + u[1] ** 2, 2 * u[0] * u[1]], [1.0, 0.5], 2),
    6: (lambda t, u: [u[0] ** 2 + math.exp(t)], [0.0], 3),
    10: (lambda t, u: [u[1], math.exp(u[0])], [0.0, 0.0], 4),
}


def neumann_1d(n, h):
    """Second-difference matrix on n nodes with ghost-point Neumann walls."""
    d = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n)).tolil()
    d[0, 1] = 2.0
    d[n - 1, n - 2] = 2.0
    return d.tocsr() / h ** 2


def pde7(n=101):
    L = neumann_1d(n, 1.0 / (n - 1))
    s = solve_ivp(lambda t, u: L @ u + u + u * u, [0, 2], np.ones(n), method="BDF",
                  jac=lambda t, u: L + sp.diags(1 + 2 * u), rtol=1e-10, atol=1e-10)
    return s.t[-1] if s.status == -1 else "no blow-up"


def pde8(n=31):
    d = neumann_1d(n, 3.0 / (n - 1))
    eye = sp.identity(n)
    L = (sp.kron(d, eye) + sp.kron(eye, d)).tocsr()
    s = solve_ivp(lambda t, u: L @ u + np.exp(u), [0, 2], np.zeros(n * n), method="BDF",
                  jac=lambda t, u: L + sp.diags(np.exp(u)), rtol=1e-10, atol=1e-10)
    return s.t[-1] if s.status == -1 else "no blow-up"


def cross9(n=199):
    h = 1.0 / (n + 1)
    x = np.arange(1, n + 1) * h
    L = (sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n)) / h ** 2).tocsr()
    s = solve_ivp(lambda t, u: L @ u, [0, 1], 4 * x * (1 - x), method="BDF", jac=L,
                  rtol=1e-10, atol=1e-10, dense_output=True)
    c = n // 2                                   # x = 1/2
    return brentq(lambda t: s.sol(t)[c] - 0.5, 0.01, 0.2, xtol=1e-15)


for k in range(1, 11):
    if k in ODE:
        f, u0, tmax = ODE[k]
        vrun("hd15.%d machine" % k, lambda: tinf(f, u0, tmax), 10, cap=60)
    else:
        vrun("hd15.%d machine" % k, {7: pde7, 8: pde8, 9: cross9}[k], 10, cap=60)
    skip("hd15.%d wp50" % k, "NDSolve (blow-up / PDE)" if k != 9 else "NDSolve (PDE)")
