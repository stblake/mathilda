"""Holdout B, family C5 (experiment 230) -- reference values.

Each problem is the DISCRETE minimax value
    E = min_p max_k |f(x_k) - p(x_k)|
(real: M equispaced points on [a, b], real p of degree n; complex: the M-th
roots of unity, complex p of degree n).  The best approximation is unique
(Haar condition), and for the complex problems -- f(conj z) = conj f(z) on a
conjugation-symmetric sample set -- it therefore has REAL coefficients, so
the references solve the real-coefficient problem.

Route A (mpmath, dps 90, this file):
  real    : the discrete Remez single-point exchange algorithm, run to
            max|e| <= |E| (1 + 1e-80);
  complex : scipy SLSQP on the epigraph form -> active set -> mpmath.findroot
            on |e_k|^2 = E^2 over the active points; certified by positive
            KKT multipliers and max(inactive)/E < 1 (convex problem => global).
Route B (independent):
  real    : scipy linprog (HiGHS simplex) on the LP form -> its alternation set
            -> PARI/gp solves the levelled (n+2)x(n+2) system at
            realprecision 120 and checks max|e| over ALL M points (bc5_ref.gp);
  complex : PARI/gp Newton (derivnum Jacobian) on the same active-set
            equations from the SLSQP point, then the all-points check.
Run:  python3 bc5_ref.py   (writes bc5_data.gp)  ;  gp -q bc5_ref.gp
"""
import os
import numpy as np
import mpmath as mp
from scipy.optimize import linprog, minimize
from scipy.special import gamma as sgamma, j0 as sj0

mp.mp.dps = 90
HERE = os.path.dirname(os.path.abspath(__file__))

# (name, f_mp, f_np, gp expression in x, a, b, n, M)
REAL = [
    ("exp", mp.exp, np.exp, "exp(x)", -1, 1, 2, 101),
    ("abs", abs, np.abs, "abs(x)", -1, 1, 4, 201),
    ("sqrt", mp.sqrt, np.sqrt, "sqrt(x)", 0, 1, 3, 101),
    ("runge", lambda x: 1 / (1 + 25 * x * x), lambda x: 1 / (1 + 25 * x * x), "1/(1+25*x^2)", -1, 1, 6, 201),
    ("gamma", mp.gamma, sgamma, "gamma(x)", 1, 2, 3, 101),
    ("atan", mp.atan, np.arctan, "atan(x)", -1, 1, 5, 201),
    ("j0", lambda x: mp.besselj(0, 10 * x), lambda x: sj0(10 * x), "besselj(0,10*x)", 0, 1, 8, 201),
]
CPLX = [
    ("exp", mp.exp, np.exp, "exp(z)", 3, 64),
    ("sqrt2", lambda z: mp.sqrt(z + 2), lambda z: np.sqrt(z + 2), "sqrt(z+2)", 2, 96),
    ("gamma2", lambda z: mp.gamma(z + 2), lambda z: sgamma(z + 2), "gamma(z+2)", 4, 128),
]


def remez(f, a, b, n, M):
    xs = [mp.mpf(a) + (mp.mpf(b) - a) * k / (M - 1) for k in range(M)]
    fs = [f(x) for x in xs]
    R = sorted(set(int(round((M - 1) * (1 - mp.cos(mp.pi * i / (n + 1))) / 2)) for i in range(n + 2)))
    assert len(R) == n + 2
    for it in range(500):
        A = mp.matrix(n + 2, n + 2)
        for i, r in enumerate(R):
            for j in range(n + 1):
                A[i, j] = xs[r] ** j
            A[i, n + 1] = (-1) ** i
        sol = mp.lu_solve(A, mp.matrix([fs[r] for r in R]))
        c, E = [sol[j] for j in range(n + 1)], sol[n + 1]
        e = [fs[k] - sum(c[j] * xs[k] ** j for j in range(n + 1)) for k in range(M)]
        m = max(range(M), key=lambda k: abs(e[k]))
        if abs(e[m]) <= abs(E) * (1 + mp.mpf(10) ** -80):
            return abs(E), it
        sg = lambda k: 1 if e[k] > 0 else -1
        if m in R:
            raise RuntimeError("exchange stalled")
        if m < R[0]:
            R = [m] + R[1:] if sg(m) == sg(R[0]) else [m] + R[:-1]
        elif m > R[-1]:
            R = R[:-1] + [m] if sg(m) == sg(R[-1]) else R[1:] + [m]
        else:
            i = max(i for i in range(n + 2) if R[i] < m)
            if sg(m) == sg(R[i]):
                R[i] = m
            else:
                R[i + 1] = m
    raise RuntimeError("no convergence")


def lp_alternation(fn, a, b, n, M):
    x = a + (b - a) * np.arange(M) / (M - 1)
    V = np.vander(x, n + 1, increasing=True)
    fx = fn(x)
    # min t  s.t.  f - V c <= t,  -(f - V c) <= t
    cost = np.zeros(n + 2); cost[-1] = 1
    A = np.block([[-V, -np.ones((M, 1))], [V, -np.ones((M, 1))]])
    bvec = np.concatenate([-fx, fx])
    r = linprog(cost, A_ub=A, b_ub=bvec, bounds=[(None, None)] * (n + 1) + [(0, None)], method="highs-ds",
                options={"primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10})
    c, t = r.x[:-1], r.x[-1]
    res = fx - V @ c
    # the error curve's extrema: one point (max |r|) per run of equal sign over
    # ALL M points, then the window of n+2 consecutive extrema whose smallest
    # |r| is largest (the LP optimum equioscillates on it)
    ext, cur = [], None
    for k in range(M):
        sk = np.sign(res[k])
        if sk == 0:
            continue
        if cur is None or sk != np.sign(res[cur[-1]]):
            ext.append(k); cur = [k]
        elif abs(res[k]) > abs(res[ext[-1]]):
            ext[-1] = k; cur.append(k)
    assert len(ext) >= n + 2, ext
    w = max(range(len(ext) - n - 1), key=lambda i: min(abs(res[k]) for k in ext[i:i + n + 2]))
    alt = ext[w:w + n + 2]
    assert min(abs(res[k]) for k in alt) > 0.99 * t, ("weak alternation", [res[k] / t for k in alt])
    return t, alt, [int(np.sign(res[k])) for k in alt]


def complex_a(f, fn, n, M):
    th = 2 * np.pi * np.arange(M) / M
    z = np.exp(1j * th); fz = fn(z)
    V = np.vander(z, n + 1, increasing=True)
    d = n + 1                                  # real coefficients (see the docstring)
    err = lambda c: fz - V @ c[:d]
    x0 = np.zeros(d + 1); x0[:d] = np.linalg.lstsq(V, fz, rcond=None)[0].real; x0[d] = 1
    cons = {"type": "ineq", "fun": lambda x: x[d] - np.abs(err(x))}
    r = minimize(lambda x: x[d], x0, constraints=[cons], method="SLSQP",
                 options={"ftol": 1e-16, "maxiter": 5000})
    E0 = r.x[d]; e = np.abs(err(r.x))
    act = sorted(set(min(k, (M - k) % M) for k in range(M) if e[k] > E0 * (1 - 3e-8)))
    zs = [mp.expjpi(mp.mpf(2 * k) / M) for k in range(M)]
    fs = [f(w) for w in zs]
    def e2(cf, k):
        return abs(fs[k] - sum(cf[j] * zs[k] ** j for j in range(n + 1))) ** 2
    assert len(act) == n + 2, ("active set size", act)
    F = lambda *u: [e2(u[:n + 1], k) - u[n + 1] ** 2 for k in act]
    sol = mp.findroot(F, [mp.mpf(v) for v in r.x[:n + 1]] + [mp.mpf(E0)])
    cf, E = [sol[j] for j in range(n + 1)], sol[n + 1]
    mx = max(mp.sqrt(e2(cf, k)) for k in range(M) if min(k, M - k) not in act)
    J = mp.matrix(n + 2, n + 2)
    h = mp.mpf(10) ** -35
    for col, k in enumerate(act):
        for c in range(n + 1):
            up = list(cf); dn = list(cf); up[c] += h; dn[c] -= h
            J[c, col] = (e2(up, k) - e2(dn, k)) / (2 * h)
        J[n + 1, col] = 1
    lam = mp.lu_solve(J, mp.matrix([0] * (n + 1) + [1]))
    return abs(E), act, list(r.x[:n + 1]) + [E0], mx / abs(E), [lam[i] for i in range(n + 2)]


if __name__ == "__main__":
    gp = ["\\\\ generated by bc5_ref.py: LP alternation sets (real) and SLSQP active sets + start points (complex)",
          "RJ = vector(%d); CJ = vector(%d);" % (len(REAL), len(CPLX))]
    for k, (nm, f, fn, gx, a, b, n, M) in enumerate(REAL):
        E, it = remez(f, a, b, n, M)
        t, alt, sg = lp_alternation(fn, a, b, n, M)
        print("bc5.%d A remez (%d exchanges) E = %s" % (k + 1, it, mp.nstr(E, 75)), flush=True)
        print("        B-lp t = %.16g, alternation %s" % (t, alt), flush=True)
        gp.append("RJ[%d] = [x -> %s, %d, %d, %d, %d, %s, %s];" % (k + 1, gx, a, b, n, M, alt, sg))
    for k, (nm, f, fn, gz, n, M) in enumerate(CPLX):
        E, act, x0, ratio, lam = complex_a(f, fn, n, M)
        kk = len(REAL) + k + 1
        print("bc5.%d A newton E = %s\n        max inactive/E = %s, KKT multipliers %s" %
              (kk, mp.nstr(E, 75), mp.nstr(ratio, 8), [mp.nstr(l, 5) for l in lam]), flush=True)
        gp.append("CJ[%d] = [z -> %s, %d, %d, %s, [%s]];" % (k + 1, gz, n, M, act,
                                                        ",".join("%.17g" % v for v in x0)))
    open(os.path.join(HERE, "bc5_data.gp"), "w").write("\n".join(gp) + "\n")
