"""Experiment 219 -- SIAM 100-Digit Challenge, holdout A, family d14 (Python column).

The 1-D two-point problems (hd14.1-3):
  machine: scipy.integrate.solve_bvp (tol 1e-10, 11-node initial mesh, zero /
           sin(pi x)/2 initial guess) inside brentq on the bracket of the .m file;
  wp50:    mpmath has no BVP solver; the idiomatic statement is shooting:
           odefun (Taylor) for the IVP from x = 0 with slope c, findroot
           (secant) on c for the right-hand boundary condition, and findroot
           (Anderson, bracketing) on the parameter over the same bracket.
The 2-D problems (hd14.4-10) are finite-element solves; neither scipy nor
mpmath has a PDE/FEM solver, so those rows are SKIP rather than a hand-written
discretisation (as in experiment 189).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, skip, require
import numpy as np
import mpmath
from scipy.integrate import solve_bvp
from scipy.optimize import brentq

require(["scipy.integrate:solve_bvp", "scipy.optimize:brentq", "mpmath:odefun", "mpmath:findroot"])

X = np.linspace(0, 1, 11)


def bvp(f, bc, guess):
    s = solve_bvp(f, bc, X, guess, tol=1e-10, max_nodes=100000)
    if s.status != 0:
        raise RuntimeError(s.message)
    return s.sol


def m1(a):
    s = bvp(lambda x, u: np.vstack([u[1], -np.exp(a * x)]), lambda l, r: np.array([l[0], r[1]]), np.zeros((2, 11)))
    return s(1.0)[0] - np.exp(a) / 4


def m2(a):
    s = bvp(lambda x, u: np.vstack([u[1], a * u[0] - 1]), lambda l, r: np.array([l[0], r[1]]), np.zeros((2, 11)))
    return s(1.0)[0] - 0.25


def m3(l):
    g = np.vstack([np.sin(np.pi * X) / 2, np.pi * np.cos(np.pi * X) / 2])
    s = bvp(lambda x, u: np.vstack([u[1], -l * np.exp(u[0])]), lambda a, b: np.array([a[0], b[0]]), g)
    return s(0.5)[0] - 0.5


def shoot(rhs, c0, right):
    """slope c with right(u(1)) = 0 for the IVP u(0) = 0, u'(0) = c; returns the odefun solution."""
    sol = lambda c: mpmath.odefun(rhs, 0, [mpmath.mpf(0), c])
    c = mpmath.findroot(lambda c: right(sol(c)(1)), (mpmath.mpf(c0), mpmath.mpf(c0) + mpmath.mpf("0.1")))
    return sol(c)


def w1(a):
    s = shoot(lambda x, u: [u[1], -mpmath.exp(a * x)], "1", lambda u: u[1])
    return s(1)[0] - mpmath.exp(a) / 4


def w2(a):
    s = shoot(lambda x, u: [u[1], a * u[0] - 1], "0.5", lambda u: u[1])
    return s(1)[0] - mpmath.mpf(1) / 4


def w3(l):
    s = shoot(lambda x, u: [u[1], -l * mpmath.exp(u[0])], "1", lambda u: u[0])
    return s(mpmath.mpf(1) / 2)[0] - mpmath.mpf(1) / 2


B = {1: (0.5, 5), 2: (1, 6), 3: (1, 3.3)}
for k, (mf, wf) in {1: (m1, w1), 2: (m2, w2), 3: (m3, w3)}.items():
    lo, hi = B[k]
    vrun("hd14.%d machine" % k, lambda: brentq(mf, lo, hi, xtol=1e-15), 10, cap=60)
    vrun("hd14.%d wp50" % k, lambda: mpmath.findroot(wf, (mpmath.mpf(lo), mpmath.mpf(hi)), solver="anderson"),
         48, cap=90, dps=60)
for k in range(4, 11):
    skip("hd14.%d machine" % k, "NDSolveValue (FEM)")
    skip("hd14.%d wp50" % k, "NDSolveValue (FEM)")
