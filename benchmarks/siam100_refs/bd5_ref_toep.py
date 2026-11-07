"""Holdout B, family D5 reference, method A for the NONNORMAL TOEPLITZ rows
(bd5.3, 5.6, 5.7, 5.8, 5.10): min over monic p of degree k of ||p(A)||_2.

Stage 1 (numpy, double): minimise the smooth convex log-sum-exp
  s0 (1 + log(sum_i exp(q (s_i/s0 - 1)))/q) of the singular values of p(A),
  q = 1e2 .. 1e7 (BFGS, analytic gradient) -- as siam100_refs/d5_ref.py.  This
  locates the minimiser to ~7 digits and shows which singular values merge.
Stage 2 (mpmath, dps 100): with J the exchange matrix, M(c) = p(A) J is
  SYMMETRIC (A Toeplitz => persymmetric) and its eigenvalues are +- the
  singular values of p(A).  The minimiser is either
    'pair':   lmax(M) = -lmin(M) (double top singular value from two distinct
              eigenvalues of M), KKT system  lmax + lmin = 0,
              w grad lmax + (1-w) grad(-lmin) = 0, 0 < w < 1  (k+1 unknowns),
    'smooth': one simple extreme eigenvalue, grad = 0  (k unknowns);
  which one is read off stage 1.  Eigenpairs from mpmath.eigsy, gradients by
  Hellmann-Feynman (q^T A^j J q), Jacobian by forward differences (h = 1e-45),
  Newton to |F| < 1e-85.
Method B: bd5_ref_toep.gp (PARI qfjacobi, central differences, 110 digits,
started from stage 1 written to bd5_<p>.start).
Usage: bd5_ref_toep.py p"""
import sys, time
import numpy as np, scipy.optimize as so, mpmath as mp

p = sys.argv[1]


def toep_np(n, d):
    return sum(v*np.eye(n, k=k) for k, v in d.items())


GR = {-1: -1, 0: 1, 1: 1, 2: 1, 3: 1}
CASES = {   # name: (n, diagonals (key = j - i), degree)
    '3': (48, GR, 4),
    '6': (40, {-1: 1, 1: mp.mpf(1)/4}, 3),
    '7': (48, GR, 1),
    '8': (24, GR, 6),
    '10': (32, {-1: 1, 1: 1, 2: 1}, 3),
}
n, dg, k = CASES[p]
A = toep_np(n, {kk: float(v) for kk, v in dg.items()})
P = [np.linalg.matrix_power(A, j) for j in range(k + 1)]
pA = lambda c: P[k] + sum(c[j]*P[j] for j in range(k))


def lse(c, q):
    U, s, Vt = np.linalg.svd(pA(c)); s0 = s[0]
    w = np.exp(q*(s/s0 - 1)); W = w.sum()
    g = np.array([np.sum((w/W)*np.einsum('ij,jk,ki->i', U.T, P[j], Vt.T)) for j in range(k)])
    return s0*(1 + np.log(W)/q), g


c = np.zeros(k)
for q in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7]:
    c = so.minimize(lse, c, args=(q,), jac=True, method='BFGS',
                    options={'gtol': 1e-12, 'maxiter': 20000}).x
Jn = np.fliplr(np.eye(n))
ev = np.linalg.eigvalsh((pA(c) @ Jn + (pA(c) @ Jn).T)/2)
print("stage 1: c =", list(c))
print("  eig(p(A)J) lowest:", ev[:3], " highest:", ev[-3:])
kind = 'pair' if abs(ev[-1] + ev[0]) < 1e-5*ev[-1] else 'smooth'
top = 'max' if ev[-1] > -ev[0] else 'min'
print("  structure:", kind, "" if kind == 'pair' else top, flush=True)
np.savetxt("bd5_%s.start" % p, c, fmt="%.17g")

mp.mp.dps = 100
Am = mp.matrix(n, n)
for i in range(n):
    for j in range(n):
        Am[i, j] = mp.mpf(dg.get(j - i, 0))
Jm = mp.matrix(n, n)
for i in range(n):
    Jm[i, n-1-i] = 1
Mj = []; Pk = mp.eye(n)
for j in range(k + 1):
    Mj.append(Pk*Jm); Pk = Pk*Am


def ext(cc):
    M = Mj[k] + sum((cc[j]*Mj[j] for j in range(k)), mp.zeros(n, n))
    M = (M + M.T)/2
    E, Q = mp.eigsy(M)
    e = [E[i] for i in range(n)]
    imax = max(range(n), key=lambda i: e[i]); imin = min(range(n), key=lambda i: e[i])
    return e[imax], Q[:, imax], e[imin], Q[:, imin]


def F(x):
    lM, qM, lm, qm = ext(x[:k])
    gM = [(qM.T*Mj[j]*qM)[0] for j in range(k)]
    gm = [-(qm.T*Mj[j]*qm)[0] for j in range(k)]
    if kind == 'pair':
        w = x[k]
        return [lM + lm] + [w*gM[j] + (1-w)*gm[j] for j in range(k)], max(lM, -lm), w
    g = gM if top == 'max' else gm
    return g, (lM if top == 'max' else -lm), None


x = [mp.mpf(float(v)) for v in c] + ([mp.mpf('0.5')] if kind == 'pair' else [])
m = len(x); h = mp.mpf(10)**-45
for it in range(30):
    t0 = time.time()
    f0, lam, w = F(x)
    nf = mp.sqrt(sum(v**2 for v in f0))
    print("Newton", it, "|F| =", mp.nstr(nf, 5), " ||p(A)|| =", mp.nstr(lam, 60),
          "" if w is None else " w = " + mp.nstr(w, 15), flush=True)
    if nf < mp.mpf(10)**-85:
        break
    Jac = mp.matrix(m, m)
    for kk in range(m):
        xp = list(x); xp[kk] += h
        fk = F(xp)[0]
        for i in range(m):
            Jac[i, kk] = (fk[i] - f0[i])/h
    dx = mp.lu_solve(Jac, mp.matrix(f0))
    x = [x[i] - dx[i] for i in range(m)]
    print("   step %.0fs" % (time.time() - t0), flush=True)
lM, qM, lm, qm = ext(x[:k])
print("check: lmax =", mp.nstr(lM, 30), " -lmin =", mp.nstr(-lm, 30))
print("bd5.%s %s" % (p, mp.nstr(lam, 80)))
