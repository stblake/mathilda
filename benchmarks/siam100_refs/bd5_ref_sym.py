"""Holdout B, family D5 reference, method A for the SYMMETRIC and SMOOTH rows.

Symmetric A (bd5.1, 5.4, 5.5, 5.9): p(A) is symmetric with eigenvalues p(l_i),
so ||p(A)||_2 = max_i |p(l_i)|: a discrete Chebyshev problem.  For monic p of
degree k (k free coefficients) the optimum is characterised (alternation
theorem on a finite set) by k+1 eigenvalues l_(1) < ... < l_(k+1) with
p(l_(m)) = (-1)^m h, |p| <= |h| on every eigenvalue.  The active set comes from
a double-precision LP (scipy HiGHS); the k+1 linear equations are then solved
in mpmath (dps 90) and the alternation certificate is CHECKED (signs and the
bound on every other eigenvalue), which proves optimality.
Eigenvalues: bd5.1 cos(j pi/6) exactly; bd5.4 the Kac/Clement matrix
(off-diagonal sqrt(j (n+1-j))) has eigenvalues -10, -8, ..., 10 exactly
(checked numerically); bd5.5 tridiag(-1, 2, -1), n = 20: 2 - 2 cos(j pi/21);
bd5.9 Hilbert 8x8: mpmath.eigsy at dps 90.
bd5.2 (Frobenius norm, smooth): exact normal equations in rationals
(fractions.Fraction), value = sqrt(rational) printed by mpmath.
Method B: bd5_ref.gp (PARI; Hilbert eigenvalues from polroots of the exact
characteristic polynomial; Kac eigenvalues from the exact charpoly)."""
import numpy as np, mpmath as mp
from fractions import Fraction
from scipy.optimize import linprog
mp.mp.dps = 90


def minimax(lams_mp, k, label):
    lam = np.array([float(x) for x in lams_mp])
    m = len(lam)
    # LP: min h s.t. -h <= lam^k + sum c_j lam^j <= h
    V = np.vstack([lam**j for j in range(k)]).T
    A_ub = np.vstack([np.hstack([V, -np.ones((m, 1))]), np.hstack([-V, -np.ones((m, 1))])])
    b_ub = np.concatenate([-lam**k, lam**k])
    r = linprog(np.r_[np.zeros(k), 1], A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)]*(k+1), method='highs')
    c = r.x[:k]; h = r.x[k]
    pv = lam**k + V @ c
    cand = sorted([i for i in range(m) if abs(pv[i]) > h*(1 - 1e-6)], key=lambda i: lam[i])
    from itertools import combinations
    for act in combinations(cand, k+1):
        sg = [np.sign(pv[i]) for i in act]
        if any(sg[t] == sg[t+1] for t in range(k)):
            continue
        # solve p(l_(m)) = s0 (-1)^m H, unknowns c_0..c_{k-1}, H
        s0 = sg[0]
        M = mp.matrix(k+1, k+1); rhs = mp.matrix(k+1, 1)
        for row, i in enumerate(act):
            for j in range(k):
                M[row, j] = lams_mp[i]**j
            M[row, k] = -s0 * (-1)**row
            rhs[row] = -lams_mp[i]**k
        x = mp.lu_solve(M, rhs)
        H = x[k]
        pvals = [lams_mp[i]**k + mp.fsum(x[j]*lams_mp[i]**j for j in range(k)) for i in range(m)]
        ok = H > 0 and all(abs(v) <= H*(1 + mp.mpf(10)**-80) for v in pvals)
        alt = all(mp.sign(pvals[act[t]]) == s0*(-1)**t for t in range(k+1))
        if ok and alt:
            break
    print("%s LP h = %.15g; alternation points %s; certificate bound ok: %s, alternation ok: %s"
          % (label, h, [float(lam[i]) for i in act], ok, alt))
    print(label, mp.nstr(H, 75), flush=True)


# bd5.1
minimax([mp.cos(j*mp.pi/6) for j in range(7)], 6, "bd5.1")
# bd5.4  Kac matrix n = 11
n = 11
K = np.zeros((n, n))
for j in range(1, n):
    K[j-1, j] = K[j, j-1] = np.sqrt(j*(n - j))
print("Kac eigenvalues (numpy):", np.round(np.linalg.eigvalsh(K), 12))
minimax([mp.mpf(v) for v in range(-10, 11, 2)], 4, "bd5.4")
# bd5.5
minimax([2 - 2*mp.cos(j*mp.pi/21) for j in range(1, 21)], 4, "bd5.5")
# bd5.9 Hilbert 8
Hm = mp.matrix(8, 8)
for i in range(8):
    for j in range(8):
        Hm[i, j] = mp.mpf(1)/(i + j + 1)
E, Q = mp.eigsy(Hm)
minimax(sorted([E[i] for i in range(8)]), 3, "bd5.9")
# bd5.2 Frobenius, Grcar 48, degree 8, exact
N = 48
G = [[Fraction(0)]*N for _ in range(N)]
for i in range(N):
    for j in range(N):
        if i == j or 1 <= j - i <= 3: G[i][j] = Fraction(1)
        elif j == i - 1: G[i][j] = Fraction(-1)
Ai = np.array([[int(v) for v in row] for row in G], dtype=object)
P = [np.identity(N, dtype=object)]
for j in range(8):
    P.append(P[-1].dot(Ai))
ip = lambda X, Y: Fraction(int((X*Y).sum()))
Gm = [[ip(P[i], P[j]) for j in range(8)] for i in range(8)]
b = [ip(P[i], P[8]) for i in range(8)]
# solve Gm c = -b exactly (Gauss)
Mx = [row[:] + [-b[i]] for i, row in enumerate(Gm)]
for col in range(8):
    piv = next(r for r in range(col, 8) if Mx[r][col] != 0)
    Mx[col], Mx[piv] = Mx[piv], Mx[col]
    for r in range(8):
        if r != col and Mx[r][col] != 0:
            f = Mx[r][col]/Mx[col][col]
            Mx[r] = [Mx[r][t] - f*Mx[col][t] for t in range(9)]
c = [Mx[i][8]/Mx[i][i] for i in range(8)]
val2 = ip(P[8], P[8]) + sum(b[i]*c[i] for i in range(8))
print("bd5.2 min ||p(A)||_F^2 =", val2, " c =", [str(x) for x in c][:3], "...")
print("bd5.2", mp.nstr(mp.sqrt(mp.mpf(val2.numerator)/val2.denominator), 75))
