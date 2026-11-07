# D5 reference, method B (and the start point for method A, d5_ref.wl).
# Stage 1 (numpy, double): minimise the smooth convex log-sum-exp
#   s0 (1 + log(sum_i exp(p (s_i/s0 - 1)))/p)  of the singular values s_i of
#   B(c) = A^8 + sum_{j<8} c_j A^j,  p = 1e2 .. 1e7  (BFGS, analytic gradient).
#   As p grows the minimiser converges to that of ||B||_2 and the top two
#   singular values merge (1766.3135..., 1766.3135...; third 1756.80).
# Stage 2 (mpmath, dps 100): with J the exchange matrix, M(c) = B(c) J is
#   symmetric (A is Toeplitz) with eigenvalues +-s_i; the double top singular
#   value is lmax(M) = -lmin(M).  Newton on the 9x9 KKT system
#     lmax + lmin = 0,   w grad(lmax) + (1-w) grad(-lmin) = 0,
#   eigenpairs from mpmath.eigsy, gradients by Hellmann-Feynman, Jacobian by
#   forward differences of the gradients (h = 1e-45).  Independent code and
#   eigensolver from d5_ref.wl (WL Eigensystem, central differences, 250 digits).
import numpy as np, scipy.optimize as so, mpmath as mp, time
n = 48
A = np.eye(n) - np.eye(n, k=-1) + np.eye(n, k=1) + np.eye(n, k=2) + np.eye(n, k=3)
P = [np.linalg.matrix_power(A, j) for j in range(9)]
pA = lambda c: P[8] + sum(c[j]*P[j] for j in range(8))
def lse(c, p):
    U, s, Vt = np.linalg.svd(pA(c)); s0 = s[0]
    w = np.exp(p*(s/s0 - 1)); W = w.sum()
    g = np.array([np.sum((w/W)*np.einsum('ij,jk,ki->i', U.T, P[j], Vt.T)) for j in range(8)])
    return s0*(1 + np.log(W)/p), g
c = np.zeros(8)
for p in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7]:
    c = so.minimize(lse, c, args=(p,), jac=True, method='BFGS',
                    options={'gtol': 1e-12, 'maxiter': 20000}).x
    print("p = %g  top singular values:" % p, np.linalg.svd(pA(c), compute_uv=False)[:3], flush=True)
print("stage-1 c:", list(c), flush=True)

mp.mp.dps = 100
Am = mp.matrix(n, n)
for i in range(n):
    for j in range(n):
        if i == j or 1 <= j - i <= 3: Am[i, j] = 1
        elif j == i - 1: Am[i, j] = -1
Jm = mp.matrix(n, n)
for i in range(n): Jm[i, n-1-i] = 1
Mj = []; Pk = mp.eye(n)
for j in range(9):
    Mj.append(Pk*Jm); Pk = Pk*Am
def ext(cc):
    M = Mj[8] + sum((cc[j]*Mj[j] for j in range(8)), mp.zeros(n, n))
    M = (M + M.T)/2                     # exact symmetry against rounding
    E, Q = mp.eigsy(M)
    ev = [E[i] for i in range(n)]
    imax = max(range(n), key=lambda i: ev[i]); imin = min(range(n), key=lambda i: ev[i])
    return ev[imax], Q[:, imax], ev[imin], Q[:, imin]
def F(x):
    lM, qM, lm, qm = ext(x[:8]); w = x[8]
    gM = [(qM.T*Mj[j]*qM)[0] for j in range(8)]
    gm = [-(qm.T*Mj[j]*qm)[0] for j in range(8)]
    return [lM + lm] + [w*gM[j] + (1-w)*gm[j] for j in range(8)], lM
x = [mp.mpf(float(v)) for v in c] + [mp.mpf('0.5')]
h = mp.mpf(10)**-45
for it in range(10):
    t0 = time.time()
    f0, lam = F(x)
    nf = mp.sqrt(sum(v**2 for v in f0))
    print("Newton", it, "|F| =", mp.nstr(nf, 5), " ||B|| =", mp.nstr(lam, 70), " w =", mp.nstr(x[8], 20), flush=True)
    if nf < mp.mpf(10)**-85: break
    Jac = mp.matrix(9, 9)
    for k in range(9):
        xp = list(x); xp[k] += h
        fk = F(xp)[0]
        for i in range(9): Jac[i, k] = (fk[i] - f0[i])/h
    dx = mp.lu_solve(Jac, mp.matrix(f0))
    x = [x[i] - dx[i] for i in range(9)]
    print("   step %.0fs" % (time.time() - t0), flush=True)
print("min ||p(A)||_2 =", mp.nstr(lam, 80))
