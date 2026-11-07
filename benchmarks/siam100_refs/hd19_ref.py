"""Reference A for holdout family hd19 (exp 224): finite variable-thickness
(compliance) design problems, the generalisation of exp 194's proxy.

Each problem: cells of side h (1-D: n intervals on [0,1]; 2-D: nx x ny squares),
c constant per cell, linear (1-D) / bilinear Q1 (2-D) elements, u = 0 on the
whole boundary, load f (1 or x),
    J(c) = F . U,  K(c) U = F,  F_i = h^d f(node_i)   (exact consistent load),
minimise J over lo <= c_e <= hi, sum c_e = mean * ncells.  J is convex in c,
so the minimum VALUE is unique.
Method A (this file): scipy SLSQP (analytic gradient -U_e^T K_e U_e) locates
the active set; the KKT system (free cells: U^T K_e U = mu; volume) is solved
by mpmath.findroot at 80 digits; the KKT sign conditions are checked; the
optimal design is written to hd19_cstar_<k>.txt for the independent duality
certificate (siam100_refs/hd19_ref.wl).
Usage: python3.11 hd19_ref.py [k ...]
"""
import sys
import numpy as np
import mpmath as mp
from scipy.optimize import minimize
from fractions import Fraction as Fr

KQ = [[4, -1, -2, -1], [-1, 4, -1, -2], [-2, -1, 4, -1], [-1, -2, -1, 4]]   # /6
K1 = [[1, -1], [-1, 1]]                                                     # /h

# k: (dim, nx, ny, lo, hi, mean, load)
PROBS = {
    1: (2, 4, 4, 1, 10, Fr(11, 2), "1"),
    2: (2, 4, 4, 1, 100, Fr(103, 4), "1"),
    3: (2, 6, 6, 1, 100, Fr(101, 2), "1"),
    4: (2, 5, 5, 1, 100, Fr(101, 2), "1"),
    5: (1, 8, 0, 1, 100, Fr(101, 2), "1"),
    6: (1, 16, 0, 1, 4, Fr(5, 2), "1"),
    7: (2, 4, 4, 1, 100, Fr(101, 2), "x"),
    8: (2, 8, 4, 1, 100, Fr(101, 2), "1"),
    9: (2, 4, 4, 1, 2, Fr(3, 2), "1"),
    10: (2, 6, 6, 1, 1000, Fr(1003, 4), "1"),
}


def build(k):
    dim, nx, ny, lo, hi, mean, load = PROBS[k]
    if dim == 1:
        h = Fr(1, nx)
        node = lambda i: i - 1 if 1 <= i <= nx - 1 else -1
        elems = [[node(a), node(a + 1)] for a in range(nx)]
        N = nx - 1
        xs = [h * (i + 1) for i in range(N)]
        F = [h * (1 if load == "1" else x) for x in xs]
        loc = [[Fr(v) / h for v in row] for row in K1]
    else:
        h = Fr(1, ny)                       # square cells: domain nx*h by ny*h
        def node(i, j):
            return (i - 1) * (ny - 1) + (j - 1) if 1 <= i <= nx - 1 and 1 <= j <= ny - 1 else -1
        elems = [[node(a, b), node(a + 1, b), node(a + 1, b + 1), node(a, b + 1)]
                 for a in range(nx) for b in range(ny)]
        N = (nx - 1) * (ny - 1)
        F = [None] * N
        for i in range(1, nx):
            for j in range(1, ny):
                F[node(i, j)] = h * h * (1 if load == "1" else h * i)
        loc = [[Fr(v, 6) for v in row] for row in KQ]
    # symmetry orbits of the cells (the problem is invariant under these maps;
    # a convex problem has a group-invariant optimum, and reducing to one
    # unknown per orbit removes the flat directions of the KKT system)
    if dim == 1:
        maps = [lambda a: a, lambda a: nx - 1 - a]
        cell = lambda a: a
        orb = lambda e: min(m(e) for m in maps)
    else:
        def imgs(e):
            a, b = divmod(e, ny)
            out = set()
            xs = [a, nx - 1 - a] if load == "1" else [a]
            for p in xs:
                for q in (b, ny - 1 - b):
                    out.add(p * ny + q)
                    if nx == ny:
                        out.add(q * ny + p) if load == "1" else None
            return out
        orb = lambda e: min(imgs(e))
    return dict(elems=elems, N=N, F=F, loc=loc, lo=lo, hi=hi, V=mean * len(elems), orb=orb)


def assemble(P, c, conv):
    N, K = P["N"], None
    K = [[conv(0)] * N for _ in range(N)]
    loc = [[conv(v) for v in row] for row in P["loc"]]
    for e, idx in enumerate(P["elems"]):
        for r, ir in enumerate(idx):
            if ir < 0:
                continue
            for s, js in enumerate(idx):
                if js >= 0:
                    K[ir][js] += c[e] * loc[r][s]
    return K


def energies(P, U, conv, zero):
    loc = [[conv(v) for v in row] for row in P["loc"]]
    out = []
    for idx in P["elems"]:
        ue = [U[i] if i >= 0 else zero for i in idx]
        out.append(sum(ue[r] * loc[r][s] * ue[s] for r in range(len(idx)) for s in range(len(idx))))
    return out


def solve(k):
    P = build(k)
    ne, lo, hi, V = len(P["elems"]), P["lo"], P["hi"], P["V"]
    Fn = np.array([float(v) for v in P["F"]])

    def J(c):
        U = np.linalg.solve(np.array(assemble(P, c, float)), Fn)
        return Fn @ U, U

    def g(c):
        return -np.array(energies(P, J(c)[1], float, 0.0))

    res = minimize(lambda c: J(c)[0], np.full(ne, float(V) / ne), jac=g, method="SLSQP",
                   bounds=[(lo, hi)] * ne,
                   constraints=[{"type": "eq", "fun": lambda c: c.sum() - float(V), "jac": lambda c: np.ones_like(c)}],
                   options={"ftol": 1e-16, "maxiter": 3000})
    c = res.x
    tol = 1e-5 * (hi - lo)
    L = [e for e in range(ne) if c[e] < lo + tol]
    H = [e for e in range(ne) if c[e] > hi - tol]
    free = [e for e in range(ne) if e not in L and e not in H]
    mp.mp.dps = 80
    F = mp.matrix([mp.mpf(v.numerator) / v.denominator for v in P["F"]])
    conv = lambda v: mp.mpf(v.numerator) / v.denominator if isinstance(v, Fr) else mp.mpf(v)

    orbits = sorted(set(P["orb"](e) for e in free))
    members = {o: [e for e in free if P["orb"](e) == o] for o in orbits}

    def full(xs):
        cv = [mp.mpf(lo)] * ne
        for e in H:
            cv[e] = mp.mpf(hi)
        for i, o in enumerate(orbits):
            for e in members[o]:
                cv[e] = xs[i]
        return cv

    def U_of(cv):
        return mp.lu_solve(mp.matrix(assemble(P, cv, conv)), F)

    en0 = -g(c)
    mu0 = float(np.mean([en0[e] for e in free])) if free else float(np.median(en0))
    SC = mp.mpf(1) / mu0

    def kkt(*a):
        xs, mu = list(a[:-1]), a[-1]
        cv = full(xs)
        en = energies(P, U_of(cv), conv, mp.mpf(0))
        return [en[o] * SC - mu for o in orbits] + [(sum(cv) - conv(V)) / ne]

    if free:
        x0 = [mp.mpf(float(np.mean([c[e] for e in members[o]]))) for o in orbits] + [mp.mpf(1)]
        sol = mp.findroot(kkt, x0, tol=mp.mpf(10) ** -140, maxsteps=60)
        sol = [sol[i] for i in range(len(x0))]
        resid = max(abs(v) for v in kkt(*sol))
        xs, mu = sol[:-1], sol[-1] / SC
    else:
        xs, mu, resid = [], None, 0
    cv = full(xs)
    U = U_of(cv)
    en = energies(P, U, conv, mp.mpf(0))
    Jv = mp.fdot(F, U)
    if mu is not None:
        signs = (all(en[e] >= mu for e in H), all(en[e] <= mu for e in L),
                 all(lo < x < hi for x in xs), "free-cell energy spread %s" % mp.nstr(max(abs(en[e] - mu) for e in free) / mu, 3),
                 "orbits %d" % len(orbits))
    else:
        signs = "no free cell"
    print("hd19.%d: SLSQP %r  lo/hi/free = %d/%d/%d  KKT resid %s  signs %s" %
          (k, res.fun, len(L), len(H), len(free), mp.nstr(resid, 3), signs), flush=True)
    print("hd19.%d J = %s" % (k, mp.nstr(Jv, 72)), flush=True)
    with open("hd19_cstar_%d.txt" % k, "w") as fh:
        fh.write("\n".join(mp.nstr(v, 78) for v in cv) + "\n")


for a in (sys.argv[1:] or [str(k) for k in PROBS]):
    solve(int(a))
