"""Experiment 190 -- SIAM 100-Digit Challenge, App. D Problem 15 (Python column).

scipy has no PDE solver, so the idiomatic scipy statement is the method of
lines written by hand in the usual way: the 5-point finite-difference
Laplacian on a uniform 63 x 63 interior grid (h = 3/64, scipy.sparse.kron),
u' = L u + exp(u), solve_ivp (BDF, sparse Jacobian, rtol = atol = 1e-10) with a
terminal event when the centre value reaches 30 (t_inf - t ~ e^-30 there);
BDF in fact gives up first (step size below spacing, at u ~ 24, t_inf - t ~
1e-10), and that time -- like the end of NDSolve's domain -- is the answer.
This is the same PDE as the WL row, but the semi-discretisation is chosen
here, not by a library, so its truncation error (O(h^2), ~3 digits) is part
of the scipy answer.  mpmath has no PDE or stiff solver: wp50 is SKIP.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp

require(["scipy.integrate:solve_ivp", "scipy.sparse:kron"])


def tinf(n=63):
    h = 3.0 / (n + 1)
    d2 = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n)) / h ** 2
    eye = sp.identity(n)
    lap = (sp.kron(d2, eye) + sp.kron(eye, d2)).tocsr()
    c = (n // 2) * n + n // 2                      # the centre node (1.5, 1.5)
    ev = lambda t, u: u[c] - 30.0
    ev.terminal = True
    sol = solve_ivp(lambda t, u: lap @ u + np.exp(u), [0, 3], np.zeros(n * n),
                    method="BDF", jac=lambda t, u: lap + sp.diags(np.exp(u)),
                    rtol=1e-10, atol=1e-10, events=ev)
    # blow-up: the event, or else the time where BDF gives up ("step size
    # less than spacing") -- the analogue of where NDSolve stops.
    return sol.t_events[0][0] if len(sol.t_events[0]) else sol.t[-1]


vrun("d15.machine", tinf, 6)
skip("d15.wp50", "NDSolve (PDE)")
