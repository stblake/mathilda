"""Reference A for holdout-B family D19 (exp 254): finite convex "variable
thickness" design problems  min_c J(c) = F . K(c)^-1 F,
K(c) = sum_e c_e K_e (assembled), lo <= c_e <= hi, sum_e c_e = V.

The ten problems (builders below; the same assembly as the experiment halves):
  bd19.1  1-D bar [0,1], 8 P1 elements, u(0)=0, free end (cantilever), f=1,  c in [1,10],  sum 32
  bd19.2  1-D bar, 12 elements, u(0)=u(1)=0, f=1,                           c in [10,100], sum 606
  bd19.3  unit square, 4x4 Q1 cells, u=0 on the boundary, f=1,              c in [1,10],  sum 48
  bd19.4  rectangle [0,2]x[0,1], 8x4 Q1 cells, u=0 on the boundary, f=1,    c in [1,100], sum 1616
  bd19.5  unit square, 6x6 Q1, u=0 on the boundary, unit POINT load at the centre node, c in [1,100], sum 1818
  bd19.6  unit square, 6x6 Q1, u=0 only on {y=0, x<=1/3} (heat sink), f=1, c in [1,100], sum 360
  bd19.7  4x4-node grid graph (24 resistors), node (0,0) grounded, unit current in at (1,2):
          J = effective resistance,                                           c in [1,10],  sum 60
  bd19.8  unit square, 4x4 cells each cut by its (a,b)-(a+1,b+1) diagonal into two P1 triangles
          sharing one c per cell (16 unknowns), u=0 on the boundary, f=1,    c in [1,100], sum 808
  bd19.9  L-shape: unit square on a 6x6 Q1 grid minus the cells a,b >= 3 (27 cells), u=0 on its
          boundary, f=1,                                                      c in [1,100], sum 1350
  bd19.10 unit square, 10x10 Q1, u=0 on the boundary, f=1,                  c in [1,100], sum 5050
f = 1 means the consistent load of a unit source (each element adds Int phi_i).

Method A (this file): scipy SLSQP (analytic gradient -E_e(U), E_e = U^T K_e U)
locates the active set; Newton's method in mpmath at 100 digits (analytic
Jacobian dE_e/dc_f = -2 (K_e U)^T K^-1 (K_f U)) solves the KKT system
E_e(U) = mu on the free set, sum c = V (the optimal DESIGN is often not unique -- a
direction d with Sum d_e K_e U = 0 makes the Jacobian singular -- and then as
many free c_e as the nullity are held at their SLSQP values; the VALUE is
unique by convexity); the KKT sign conditions (E_e <= mu at c = lo, >= mu at
c = hi) are checked.  Writes bd19_<k>.gpdat for method B
(bd19_ref.gp: an exact-rational duality certificate in PARI).
Usage: python3.11 bd19_ref.py [k ...]
"""
import sys
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

KQ = [[4, -1, -2, -1], [-1, 4, -1, -2], [-2, -1, 4, -1], [-1, -2, -1, 4]]   # x 1/6, corners ccw


def bar(n, clamp_right):
    h = Fr(1, n)
    dof = {i: i - 1 for i in range(1, n + (0 if clamp_right else 1))}
    N = len(dof)
    F = [Fr(0)] * N
    els = []
    for e in range(n):
        nodes = [e, e + 1]
        els.append(([dof.get(v, -1) for v in nodes], [[1 / h, -1 / h], [-1 / h, 1 / h]]))
        for v in nodes:
            if v in dof:
                F[dof[v]] += h / 2
    return N, els, F


def q1(nx, ny, h, cell_ok, node_free, point=None):
    cells = [(a, b) for a in range(nx) for b in range(ny) if cell_ok(a, b)]
    used = sorted({v for a, b in cells for v in ((a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1))})
    free = [v for v in used if node_free(*v)]
    dof = {v: i for i, v in enumerate(free)}
    F = [Fr(0)] * len(free)
    els = []
    for a, b in cells:
        nodes = [(a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1)]
        els.append(([dof.get(v, -1) for v in nodes], [[Fr(k, 6) for k in r] for r in KQ]))
        if point is None:
            for v in nodes:
                if v in dof:
                    F[dof[v]] += h * h / 4
    if point is not None:
        F[dof[point]] = Fr(1)
    return len(free), els, F


def interior(nx, ny):
    return lambda i, j: 0 < i < nx and 0 < j < ny


def p1cell(n):
    """each cell = two P1 right triangles split along its (a,b)-(a+1,b+1) diagonal, ONE c per
    cell: the summed cell matrix (1/2)[[2,-1,0,-1],[-1,2,-1,0],[0,-1,2,-1],[-1,0,-1,2]] (the
    5-point stencil) and the load h^2 (1/3, 1/6, 1/3, 1/6) on the corners (a,b),(a+1,b),
    (a+1,b+1),(a,b+1)"""
    h = Fr(1, n)
    dof = {(i, j): (i - 1) * (n - 1) + (j - 1) for i in range(1, n) for j in range(1, n)}
    F = [Fr(0)] * len(dof)
    Kc = [[Fr(x, 2) for x in r] for r in ([2, -1, 0, -1], [-1, 2, -1, 0], [0, -1, 2, -1], [-1, 0, -1, 2])]
    w = [h * h / 3, h * h / 6, h * h / 3, h * h / 6]
    els = []
    for a in range(n):
        for b in range(n):
            nodes = [(a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1)]
            els.append(([dof.get(v, -1) for v in nodes], Kc))
            for v, wv in zip(nodes, w):
                if v in dof:
                    F[dof[v]] += wv
    return len(dof), els, F


def grid_graph(m):
    nodes = [(i, j) for i in range(m) for j in range(m)]
    dof = {v: k - 1 for k, v in enumerate(nodes) if k > 0}       # node (0,0) grounded
    els = []
    for i in range(m):
        for j in range(m):
            for nb in ((i + 1, j), (i, j + 1)):
                if nb[0] < m and nb[1] < m:
                    els.append(([dof.get((i, j), -1), dof.get(nb, -1)], [[Fr(1), Fr(-1)], [Fr(-1), Fr(1)]]))
    F = [Fr(0)] * len(dof)
    F[dof[(1, 2)]] = Fr(1)                                       # unit current in at (1,2)
    return len(dof), els, F


def problem(k):
    if k == 1: return bar(8, False) + (1, 10, 32)
    if k == 2: return bar(12, True) + (10, 100, 606)
    if k == 3: return q1(4, 4, Fr(1, 4), lambda a, b: True, interior(4, 4)) + (1, 10, 48)
    if k == 4: return q1(8, 4, Fr(1, 4), lambda a, b: True, interior(8, 4)) + (1, 100, 1616)
    if k == 5: return q1(6, 6, Fr(1, 6), lambda a, b: True, interior(6, 6), point=(3, 3)) + (1, 100, 1818)
    if k == 6: return q1(6, 6, Fr(1, 6), lambda a, b: True, lambda i, j: not (j == 0 and i <= 2)) + (1, 100, 360)
    if k == 7: return grid_graph(4) + (1, 10, 60)
    if k == 8: return p1cell(4) + (1, 100, 808)
    if k == 9:
        ok = lambda a, b: not (a >= 3 and b >= 3)
        cell = lambda i, j: all(ok(a, b) and 0 <= a < 6 and 0 <= b < 6 for a in (i - 1, i) for b in (j - 1, j))
        return q1(6, 6, Fr(1, 6), ok, cell) + (1, 100, 1350)
    if k == 10: return q1(10, 10, Fr(1, 10), lambda a, b: True, interior(10, 10)) + (1, 100, 5050)


def assemble(N, els, c, zero):
    K = [[zero] * N for _ in range(N)]
    for (idx, Ke), ce in zip(els, c):
        for r, ir in enumerate(idx):
            if ir < 0: continue
            for s, js in enumerate(idx):
                if js >= 0:
                    K[ir][js] += ce * Ke[r][s]
    return K


def energies(els, U, zero):
    out = []
    for idx, Ke in els:
        ue = [U[i] if i >= 0 else zero for i in idx]
        out.append(sum(ue[r] * Ke[r][s] * ue[s] for r in range(len(idx)) for s in range(len(idx))))
    return out


def solve_A(k, verbose=True):
    N, els, Fq, lo, hi, V = problem(k)
    m = len(els)
    elsf = [(idx, [[float(x) for x in r] for r in Ke]) for idx, Ke in els]
    Ff = np.array([float(x) for x in Fq])

    def JU(c):
        K = np.array(assemble(N, elsf, c, 0.0))
        U = np.linalg.solve(K, Ff)
        return Ff @ U, U

    fun = lambda c: JU(c)[0]
    jac = lambda c: -np.array(energies(elsf, JU(c)[1], 0.0))
    res = minimize(fun, np.full(m, V / m), jac=jac, method="SLSQP", bounds=[(lo, hi)] * m,
                   constraints=[{"type": "eq", "fun": lambda c: c.sum() - V, "jac": lambda c: np.ones(m)}],
                   options={"ftol": 1e-16, "maxiter": 3000})
    c = res.x
    tol = 1e-6 * (hi - lo)
    LO = {e for e in range(m) if c[e] < lo + tol}
    # an element with no free node carries no energy (E_e = 0 < mu): c_e = lo at the optimum
    LO |= {e for e in range(m) if all(i < 0 for i in els[e][0])}
    HI = {e for e in range(m) if c[e] > hi - tol} - LO
    if verbose:
        print("bd19.%d: N=%d m=%d  SLSQP J=%r  (%s)  lo %d hi %d free %d" % (k, N, m, res.fun, res.message, len(LO), len(HI), m - len(LO) - len(HI)), flush=True)
    mp.mp.dps = 100
    Fm = mp.matrix([mp.mpf(x.numerator) / x.denominator for x in Fq])
    elsm = [(idx, [[mp.mpf(x.numerator) / x.denominator for x in r] for r in Ke]) for idx, Ke in els]
    for it_outer in range(10):
        free = [e for e in range(m) if e not in LO and e not in HI]
        cv = [mp.mpf(lo) if e in LO else mp.mpf(hi) if e in HI else mp.mpf(c[e]) for e in range(m)]
        # Non-unique optimal DESIGN: J* and the optimal state U* are unique (strict concavity
        # of the dual), but the optimal designs form a polytope {c : K(c) U* = F, bounds, sum}
        # whenever d -> (Sum_e d_e K_e U, Sum_e d_e) has a null space on the free set; the KKT
        # Jacobian is then singular.  Fix that many free c_e (the nullity) at their SLSQP values
        # -- chosen by pivoted QR of the null basis -- and solve the KKT system for the rest:
        # a point of the optimal polytope, with the same J*.  Checked: E_e = mu also holds on
        # the fixed elements (to the working precision), and method B certifies the value.
        Uf = JU(np.array([float(x) for x in cv]))[1]
        fixed = []
        if free:
            cols = []
            for e in free:
                idx, Ke = elsf[e]
                v = np.zeros(N)
                ue = np.array([Uf[i] if i >= 0 else 0.0 for i in idx])
                for r, ir in enumerate(idx):
                    if ir >= 0:
                        v[ir] = np.dot(Ke[r], ue)
                cols.append(np.concatenate([v, [np.linalg.norm(v)]]))
            Amat = np.array(cols).T
            _, sv, Vt = np.linalg.svd(Amat)
            rank = int((sv > 1e-9 * sv[0]).sum())
            null = len(free) - rank
            if null > 0:
                import scipy.linalg as _sl
                piv = _sl.qr(Vt[rank:], pivoting=True)[2]
                fixed = [free[i] for i in piv[:null]]
                if verbose:
                    print("   design non-unique: nullity %d, fixing %s" % (null, fixed), flush=True)
        free = [e for e in free if e not in fixed]
        mu = mp.mpf(float(np.mean([energies(elsf, Uf, 0.0)[e] for e in free + fixed]))) if free else mp.mpf(0)
        SC = 1 / mu if mu else mp.mpf(1)
        for it in range(60):
            K = mp.matrix(assemble(N, elsm, cv, mp.mpf(0)))
            Ki = mp.inverse(K)
            U = Ki * Fm
            E = energies(elsm, U, mp.mpf(0))
            R = [(E[e] - mu) * SC for e in free] + [sum(cv) - V]
            nr = max(abs(x) for x in R)
            if nr < mp.mpf(10) ** -90:
                break
            # global vectors K_f U
            KU = {}
            for f in free:
                idx, Ke = elsm[f]
                vec = mp.matrix(N, 1)
                ue = [U[i] if i >= 0 else 0 for i in idx]
                for r, ir in enumerate(idx):
                    if ir >= 0:
                        vec[ir] = sum(Ke[r][s] * ue[s] for s in range(len(idx)))
                KU[f] = vec
            W = {f: Ki * KU[f] for f in free}
            nf = len(free)
            Jm = mp.matrix(nf + 1, nf + 1)
            for a, e in enumerate(free):
                for b, f in enumerate(free):
                    Jm[a, b] = -2 * SC * sum(KU[e][i] * W[f][i] for i in range(N))
                Jm[a, nf] = -SC
                Jm[nf, a] = 1
            dz = mp.lu_solve(Jm, mp.matrix(R))
            for a, e in enumerate(free):
                cv[e] -= dz[a]
            mu -= dz[nf]
        bad_lo = [e for e in LO if E[e] > mu]
        bad_hi = [e for e in HI if E[e] < mu]
        bad_free = [e for e in free + fixed if cv[e] < lo or cv[e] > hi]
        if fixed and verbose:
            print('   fixed elements: max |E_e/mu - 1| = %s' % mp.nstr(max(abs(E[e] / mu - 1) for e in fixed), 3), flush=True)
        if verbose:
            print("   Newton it=%d |R|=%s  KKT violations lo %d hi %d bounds %d" % (it, mp.nstr(nr, 3), len(bad_lo), len(bad_hi), len(bad_free)), flush=True)
        if not (bad_lo or bad_hi or bad_free):
            break
        LO -= set(bad_lo); HI -= set(bad_hi)
        for e in bad_free:
            (LO if cv[e] < lo else HI).add(e)
    Jv = sum(Fm[i] * U[i] for i in range(N))
    if verbose:
        print("bd19.%d A  J* = %s" % (k, mp.nstr(Jv, 75)), flush=True)
    with open("bd19_%d.gpdat" % k, "w") as fh:
        fh.write("N = %d; m = %d; lo = %d; hi = %d; V = %d;\n" % (N, m, lo, hi, V))
        fh.write("F = [%s]~;\n" % ", ".join(str(x) for x in Fq))
        fh.write("IDX = [%s];\n" % ", ".join("[%s]" % ", ".join(str(i + 1) for i in idx) for idx, _ in els))
        fh.write("KE = [%s];\n" % ", ".join("[%s]" % "; ".join(", ".join(str(x) for x in r) for r in Ke) for _, Ke in els))
        fh.write("CS = [%s];\n" % ", ".join('"%s"' % mp.nstr(x, 88) for x in cv))
        fh.write("FREE = [%s];\n" % ", ".join(str(e + 1) for e in free))
    return Jv


if __name__ == "__main__":
    for a in (sys.argv[1:] or [str(k) for k in range(1, 11)]):
        solve_A(int(a))
