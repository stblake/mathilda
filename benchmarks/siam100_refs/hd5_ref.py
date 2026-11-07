"""Holdout A, family hd5 (experiment 210): min over a normalised polynomial p
of ||p(A)||, method B (mpmath).  Problems as in matrix_polynomial_family.m.

Write B(c) = B0 + sum_j c_j Bj (Bj = A^j, B0 = A^k for monic p; B0 = I,
j = 1..k for p(0) = 1).
  * Frobenius norm: a linear least-squares problem; exact normal equations
    (Gram matrix of the Bj) solved at dps 100.
  * A normal (diagonal or symmetric): ||B(c)||_2 = max_i |p(mu_i)| over the
    eigenvalues mu_i of A -- a discrete Chebyshev problem, solved as the max
    over (k+1)-point subsets of the alternation level (de la Vallee Poussin),
    at dps 100, checked |p| <= t at every eigenvalue.
  * A Toeplitz, nonnormal: M(c) = B(c) J is symmetric (J the exchange matrix)
    with |eig M| = singular values of B.  Stage 1: smoothing continuation of
    log-sum-exp of +-eig M (BFGS).  Stage 2: Newton at dps 80 on the KKT system
    of the active simple eigenvalues sigma_i lambda_i(c) = t,
    sum_i w_i sigma_i grad lambda_i = 0, sum w_i = 1 (Hellmann-Feynman
    gradients, forward-difference Jacobian); dual weights w_i > 0 checked.
Usage: hd5_ref.py [--stage1] k [k ...]   (--stage1: only write hd5_start.wl for hd5_ref.wl)"""
import sys, time
import numpy as np, scipy.optimize as so, mpmath as mp
from fractions import Fraction as Fr


def toeplitz(n, diags):
    """diags: {offset: value}, offset = j - i."""
    A = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        for off, v in diags.items():
            j = i + off
            if 0 <= j < n:
                A[i][j] = Fr(v)
    return A


def grcar(n, s):
    d = {-1: -1, 0: 1}
    for j in range(1, s + 1):
        d[j] = 1
    return toeplitz(n, d)


# (kind, data, k, norm, normalisation); kind 'diag' -> eigenvalue list (mp),
# 'sym' -> symmetric (eigenvalues given in closed form), 'toep' -> Fraction matrix
def eig_cheb(j):
    return mp.cos(j * mp.pi / 8)


PROBLEMS = {
    1: ('diag', lambda: [mp.cos(j * mp.pi / 8) for j in range(9)], 8, '2', 'monic'),
    2: ('diag', lambda: [mp.mpf(j) for j in range(1, 11)], 4, '2', 'monic'),
    3: ('toep', lambda: grcar(48, 3), 8, 'fro', 'monic'),
    4: ('diag', lambda: [2 - 2 * mp.cos(j * mp.pi / 21) for j in range(1, 21)], 5, '2', 'monic'),
    5: ('toep', lambda: grcar(32, 3), 6, '2', 'monic'),
    6: ('toep', lambda: grcar(40, 2), 6, '2', 'monic'),
    7: ('diag', lambda: [mp.mpf(j) for j in range(1, 11)], 3, '2', 'gmres'),
    8: ('toep', lambda: toeplitz(12, {0: 1, 1: 2}), 4, '2', 'monic'),
    9: ('toep', lambda: toeplitz(24, {-1: Fr(1, 4), 1: 1}), 6, '2', 'monic'),
    10: ('toep', lambda: grcar(64, 3), 10, '2', 'monic'),
}


def matmul(A, B):
    n = len(A)
    return [[sum(A[i][l] * B[l][j] for l in range(n) if A[i][l]) for j in range(n)] for i in range(n)]


def basis(A, k, normal):
    n = len(A)
    I = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    P = [I]
    for _ in range(k):
        P.append(matmul(P[-1], A))
    if normal == 'monic':
        return P[k], P[:k]
    return P[0], P[1:k + 1]


def solve_fro(B0, Bs):
    q = lambda x: mp.mpf(x.numerator) / x.denominator
    ip = lambda X, Y: sum(X[i][j] * Y[i][j] for i in range(len(X)) for j in range(len(X)))
    G = mp.matrix([[q(ip(a, b)) for b in Bs] for a in Bs])
    g = mp.matrix([q(ip(a, B0)) for a in Bs])
    c = mp.lu_solve(G, -g)
    val2 = q(ip(B0, B0)) + sum(c[i] * g[i] for i in range(len(Bs)))
    return mp.sqrt(val2)


def solve_diag(mus, k, normal):
    m = len(mus)
    # monomials: monic -> unknowns c_0..c_{k-1}, fixed x^k; gmres -> c_1..c_k, fixed 1
    pw = list(range(k)) if normal == 'monic' else list(range(1, k + 1))
    fixed = k if normal == 'monic' else 0
    # de la Vallee Poussin / discrete alternation: the minimax value on the m
    # points is the MAX over (k+1)-point subsets of the subset's alternation
    # level |h| (p(x_i) = (-1)^i h on the sorted subset).  Rank the subsets in
    # doubles, then solve the winner at dps 100 and check |p| <= |h| everywhere.
    from itertools import combinations
    nv = len(pw) + 1
    mf = np.array([float(x) for x in mus])
    best = []
    for S in combinations(range(m), nv):
        xs = mf[list(S)]
        Mm = np.array([[x ** j for j in pw] + [-(-1) ** r] for r, x in enumerate(xs)])
        try:
            sol = np.linalg.solve(Mm, -xs ** fixed)
        except np.linalg.LinAlgError:
            continue
        best.append((abs(sol[-1]), S))
    best.sort(reverse=True)
    for _, S in best[:20]:
        Mm = mp.matrix([[mus[i] ** j for j in pw] + [-(-1) ** r] for r, i in enumerate(S)])
        sol = mp.lu_solve(Mm, mp.matrix([-mus[i] ** fixed for i in S]))
        cc, tt = [sol[i] for i in range(len(pw))], sol[len(pw)]
        allv = [sum(cc[i] * x ** pw[i] for i in range(len(pw))) + x ** fixed for x in mus]
        if max(abs(v) for v in allv) <= abs(tt) * (1 + mp.mpf(10) ** -60):
            print("   alternation set", list(S), flush=True)
            return abs(tt)
    raise RuntimeError("no feasible alternation set")


def solve_toep(B0, Bs, dps):
    n = len(B0); k = len(Bs)
    J = np.fliplr(np.eye(n))
    Mf = [np.array(B, dtype=float) @ J for B in [B0] + Bs]
    Mof = lambda c: Mf[0] + sum(c[j] * Mf[j + 1] for j in range(k))

    def lse(c, p):
        E, Q = np.linalg.eigh((Mof(c) + Mof(c).T) / 2)
        s = np.concatenate([E, -E]); s0 = s.max()
        w = np.exp(p * (s / s0 - 1)); W = w.sum()
        g = np.zeros(k)
        for j in range(k):
            d = np.einsum('ij,jk,ki->i', Q.T, Mf[j + 1], Q)
            g[j] = np.sum(w[:n] / W * d) - np.sum(w[n:] / W * d)
        return s0 * (1 + np.log(W) / p), g
    c = np.zeros(k)
    for p in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7]:
        c = so.minimize(lse, c, args=(p,), jac=True, method='BFGS',
                        options={'gtol': 1e-13, 'maxiter': 50000}).x
    E = np.linalg.eigvalsh(Mof(c))
    top = max(abs(E))
    act = [i for i in range(n) if abs(abs(E[i]) - top) < 1e-5 * top]
    sg = [1 if E[i] > 0 else -1 for i in act]
    print("   stage 1: |eig| top", sorted(abs(E))[::-1][:5], "active", act, sg, flush=True)
    assert len(set(sg)) == len(sg), "repeated active eigenvalue of one sign"
    STAGE1.append(c)
    if STAGE1_ONLY:
        return mp.mpf(float(top))
    mp.mp.dps = dps
    Jm = mp.matrix(n, n)
    for i in range(n):
        Jm[i, n - 1 - i] = 1
    Mm = [mp.matrix([[mp.mpf(x.numerator) / x.denominator for x in row] for row in B]) * Jm
          for B in [B0] + Bs]
    r = len(act)

    def eigs(cc):
        M = Mm[0] + sum((cc[j] * Mm[j + 1] for j in range(k)), mp.zeros(n, n))
        M = (M + M.T) / 2
        Ev, Q = mp.eigsy(M)
        order = sorted(range(n), key=lambda i: Ev[i])
        return [Ev[i] for i in order], [Q[:, i] for i in order]

    # active eigenvalues are identified by their rank in the sorted spectrum
    ranks = [sorted(range(n), key=lambda i: E[i]).index(a) for a in act]

    def F(x):
        cc, t, w = x[:k], x[k], x[k + 1:]
        Ev, Qv = eigs(cc)
        out = [sg[a] * Ev[ranks[a]] - t for a in range(r)]
        for j in range(k):
            out.append(sum(w[a] * sg[a] * (Qv[ranks[a]].T * Mm[j + 1] * Qv[ranks[a]])[0]
                           for a in range(r)))
        out.append(sum(w) - 1)
        return out
    x = [mp.mpf(float(v)) for v in c] + [mp.mpf(float(top))] + [mp.mpf(1) / r] * r
    h = mp.mpf(10) ** (-dps // 2 + 5)
    for it in range(25):
        f0 = F(x)
        nf = mp.sqrt(sum(v ** 2 for v in f0))
        print("   Newton", it, "|F| =", mp.nstr(nf, 5), "t =", mp.nstr(x[k], 30), flush=True)
        if nf < mp.mpf(10) ** (-dps + 15):
            break
        m = len(x)
        Jac = mp.matrix(m, m)
        for q in range(m):
            xp = list(x); xp[q] += h
            fq = F(xp)
            for i in range(m):
                Jac[i, q] = (fq[i] - f0[i]) / h
        dx = mp.lu_solve(Jac, mp.matrix(f0))
        x = [x[i] - dx[i] for i in range(m)]
    print("   dual weights", [mp.nstr(v, 10) for v in x[k + 1:]], flush=True)
    assert all(v > 0 for v in x[k + 1:])
    Ev, _ = eigs(x[:k])
    assert max(abs(v) for v in Ev) <= x[k] * (1 + mp.mpf(10) ** -50)
    return x[k]


STAGE1, STAGE1_ONLY = [], False

if __name__ == "__main__":
    if sys.argv[1:2] == ["--stage1"]:
        STAGE1_ONLY = True; sys.argv.pop(1)
        fh = open("hd5_start.wl", "w")
    for kk in [int(a) for a in sys.argv[1:]] or sorted(PROBLEMS):
        kind, mk, k, nrm, normal = PROBLEMS[kk]
        t0 = time.time()
        mp.mp.dps = 100
        if kind == 'diag':
            v = solve_diag(mk(), k, normal)
        elif nrm == 'fro':
            v = solve_fro(*basis(mk(), k, normal))
        else:
            v = solve_toep(*basis(mk(), k, normal), dps=80)
        print("RESULT", kk, mp.nstr(v, 70), "(%.0fs)" % (time.time() - t0), flush=True)
        if STAGE1_ONLY and STAGE1:
            fh.write("hd5start[%d] = {%s};\n" % (kk, ", ".join(repr(float(x)).replace('e', '*^') for x in STAGE1.pop())))
    if STAGE1_ONLY:
        fh.close()
