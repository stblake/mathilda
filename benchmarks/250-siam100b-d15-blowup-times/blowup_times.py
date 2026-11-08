"""Experiment 250 -- SIAM 100-Digit Challenge holdout B for App. D Problem 15 (Python column).

scipy has no PDE solver, so (as in exp. 190) the idiomatic scipy statement of
the PDE rows is the method of lines written by hand in the usual way: 2nd-order
finite differences (cell-vertex grid; Neumann ends by reflection; the radial
Laplacian u_rr + u_r/r with 2 u_rr at r = 0 for the radially symmetric disk
problem bd15.9), solve_ivp BDF with the sparsity pattern of the Jacobian,
rtol = atol = 1e-10, and a terminal event when max u passes a threshold that
lies within ~1e-12 of the singularity (u = 30 for the e^u / log-type
blow-ups, 1e12 for the pole-type ones, 1 - 1e-8 for the quenching row);
otherwise the time at which BDF gives up -- the analogue of where NDSolve
stops.  The semi-discretisation is chosen here, not by a library, so its
truncation error is part of the scipy answer (the constant-data Neumann rows
1, 5, 8 are exact on any grid).  The ODE rows use DOP853 (rtol = atol = 1e-13)
with the same kind of threshold event.  mpmath has no PDE / stiff solver and
no singularity detection: every wp50 row is SKIP.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp

require(["scipy.integrate:solve_ivp", "scipy.sparse:kron"])


def d2_neumann(n, h):
    """(n+1)x(n+1) second difference with reflecting (zero-flux) ends."""
    m = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n + 1, n + 1)).tolil()
    m[0, 1] = 2.0
    m[n, n - 1] = 2.0
    return m.tocsr() / h ** 2


def d1_neumann(n, h):
    m = sp.diags([-1.0, 1.0], [-1, 1], shape=(n + 1, n + 1)).tolil()
    m[0, 1] = 0.0
    m[n, n - 1] = 0.0
    return m.tocsr() / (2 * h)


def d2_dirichlet(n, h):
    return sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n)).tocsr() / h ** 2


def d1_dirichlet(n, h):
    return sp.diags([-1.0, 1.0], [-1, 1], shape=(n, n)).tocsr() / (2 * h)


def run(rhs, y0, tmax, top, pattern, method="BDF", tol=1e-10):
    ev = lambda t, u: np.max(u) - top
    ev.terminal = True
    kw = dict(jac_sparsity=pattern) if method == "BDF" else {}
    s = solve_ivp(rhs, [0, tmax], y0, method=method, rtol=tol, atol=tol, events=ev, **kw)
    return s.t_events[0][0] if len(s.t_events[0]) else s.t[-1]


def pat(*ms):
    p = ms[0].copy()
    for m in ms[1:]:
        p = p + m
    return (abs(p) + sp.identity(p.shape[0])) != 0


def neumann_1d(f, n, a, b, u0, tmax, top, grad=None):
    h = (b - a) / n
    D2, D1 = d2_neumann(n, h), d1_neumann(n, h)
    x = np.linspace(a, b, n + 1)
    if grad is None:
        rhs = lambda t, u: D2 @ u + f(u)
    else:
        rhs = lambda t, u: D2 @ u + grad(u, D1 @ u)
    return run(rhs, u0(x), tmax, top, pat(D2))


def case1():
    return neumann_1d(lambda u: u * u + u, 50, 0.0, 1.0, lambda x: np.ones_like(x), 2.0, 1e12)


def case2():
    n = 63; h = np.pi / (n + 1)
    D2, D1 = d2_dirichlet(n, h), d1_dirichlet(n, h)
    rhs = lambda t, u: D2 @ u - (D1 @ u) ** 2 + np.exp(u)
    return run(rhs, np.zeros(n), 3.0, 30.0, pat(D2))


def case3():
    n = 40; h = np.pi / n
    D2, D1 = d2_neumann(n, h), d1_neumann(n, h); I = sp.identity(n + 1)
    L = sp.kron(D2, I) + sp.kron(I, D2); Gx = sp.kron(D1, I); Gy = sp.kron(I, D1)
    x = np.linspace(0, np.pi, n + 1); X, Y = np.meshgrid(x, x, indexing="ij")
    rhs = lambda t, u: L @ u - (Gx @ u) ** 2 - (Gy @ u) ** 2 + np.exp(u)
    return run(rhs, (-np.log(2 + np.cos(X) * np.cos(Y))).ravel(), 3.0, 30.0, pat(L, Gx, Gy))


def case4():
    return neumann_1d(None, 64, 0.0, np.pi, lambda x: 1 / (2 + np.cos(x)), 3.0, 1e12,
                      grad=lambda u, ux: -2 * ux * ux / u + u * u)


def case5():
    return neumann_1d(lambda u: np.exp(u) / (1 - u), 50, 0.0, 1.0, lambda x: np.zeros_like(x), 1.0, 1 - 1e-8)


def case6():
    return run(lambda t, y: [y[0] ** 2 + t * t], [0.0], 3.0, 1e15, None, method="DOP853", tol=1e-13)


def case7():
    ev = lambda t, y: y[0] - 1e30
    ev.terminal = True
    s = solve_ivp(lambda t, y: [y[1], 6 * y[0] ** 2], [0, 2], [1.0, 0.0], method="DOP853",
                  rtol=1e-13, atol=1e-13, events=ev)
    return s.t_events[0][0] if len(s.t_events[0]) else s.t[-1]


def case8():
    n = 50; h = 1.0 / n
    D2 = d2_neumann(n, h); Z = sp.csr_matrix((n + 1, n + 1))
    rhs = lambda t, w: np.concatenate([D2 @ w[:n + 1] + w[n + 1:] ** 2, D2 @ w[n + 1:] + w[:n + 1] ** 2])
    P = sp.bmat([[D2, sp.identity(n + 1)], [sp.identity(n + 1), D2]])
    y0 = np.concatenate([np.ones(n + 1), 0.5 * np.ones(n + 1)])
    return run(rhs, y0, 2.0, 1e12, pat(P))


def case9():
    n = 128; R = 3.0; h = R / n
    r = np.arange(n) * h                      # nodes 0..n-1; u(R) = 0 at node n
    main = -2.0 * np.ones(n); lo = np.ones(n - 1); up = np.ones(n - 1)
    D2 = sp.diags([lo, main, up], [-1, 0, 1]).tolil(); D2[0, 1] = 2.0
    D2 = D2.tocsr() / h ** 2
    D1 = sp.diags([-np.ones(n - 1), np.ones(n - 1)], [-1, 1]).tolil(); D1[0, 1] = 0.0
    D1 = D1.tocsr() / (2 * h)
    inv_r = np.zeros(n); inv_r[1:] = 1 / r[1:]
    def rhs(t, u):
        ur = D1 @ u
        lap = D2 @ u + inv_r * ur
        lap[0] = 2 * (D2 @ u)[0]
        return lap - ur * ur + np.exp(u)
    return run(rhs, np.zeros(n), 2.0, 30.0, pat(D2, D1))


def case10():
    n = 47; h = 4.0 / (n + 1)
    D2, D1 = d2_dirichlet(n, h), d1_dirichlet(n, h); I = sp.identity(n)
    L = sp.kron(D2, I) + sp.kron(I, D2); Gx = sp.kron(D1, I); Gy = sp.kron(I, D1)
    rhs = lambda t, u: L @ u - (Gx @ u) ** 2 - (Gy @ u) ** 2 + np.exp(u)
    return run(rhs, np.zeros(n * n), 3.0, 30.0, pat(L, Gx, Gy))


CASES = [case1, case2, case3, case4, case5, case6, case7, case8, case9, case10]
for k, fn in enumerate(CASES, 1):
    vrun("bd15.%d machine" % k, fn, 10)
for k in range(1, 11):
    skip("bd15.%d wp50" % k, "NDSolve (singularity detection, MPFR)")
