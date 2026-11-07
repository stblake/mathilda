"""C7 reference: (A^-1)_11, A = diag(primes) + 1 at |i-j| = 1,2,4,...,2^k < n.
Method 1 (this script, `refine`): mixed-precision iterative refinement.  The
correction is solved in doubles (scipy cg, Jacobi-preconditioned;
splu fill-in at n = 20000 is prohibitive); the residual
b - A x is computed EXACTLY in Python integers (A is an integer matrix; x is
held as an integer scaled by 2^P), so each sweep gains ~15 digits.
Method 2 (`cg`): Jacobi-preconditioned conjugate gradients carried out in
exact fixed-point integer arithmetic with scale 2^P (no double anywhere).
Usage: c7_ref.py n method"""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from fractions import Fraction
from mpmath import mp, mpf, nstr
n = int(sys.argv[1]); method = sys.argv[2]
def primes(n):
    N = 300000 if n > 2000 else 20000
    s = np.ones(N, bool); s[:2] = False
    for i in range(2, int(N**.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0][:n]]
P = primes(n)
offs = []; k = 1
while k < n: offs.append(k); k *= 2
nbr = [[j for o in offs for j in (i-o, i+o) if 0 <= j < n] for i in range(n)]
def matvec(x):   # exact integer A x
    return [P[i]*x[i] + sum(x[j] for j in nbr[i]) for i in range(n)]
PB = 300                  # fixed-point scale 2^PB (~90 digits)
S = 1 << PB
mp.dps = 80
if method == "refine":
    rows, cols, vals = [], [], []
    for i in range(n):
        rows.append(i); cols.append(i); vals.append(float(P[i]))
        for j in nbr[i]: rows.append(i); cols.append(j); vals.append(1.0)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    Minv = sp.diags([1.0/p for p in P])   # Jacobi preconditioner
    solve = lambda r: sla.cg(A, r, rtol=1e-15, atol=0.0, M=Minv, maxiter=1000)[0]
    X = [0]*n                 # x = X / 2^PB
    for it in range(8):
        AX = matvec(X)
        R = [(S if i == 0 else 0) - AX[i] for i in range(n)]   # residual * 2^PB, exact
        m = max(abs(r) for r in R)
        if m == 0: break
        e = m.bit_length() - 60
        rf = np.array([float(Fraction(r, 1 << e)) if e > 0 else float(r) for r in R])
        d = solve(rf)
        X = [X[i] + (int(round(d[i] * 2**40)) << e) // (1 << 40) if e > 40 else X[i] + int(round(d[i] * 2.0**e)) for i in range(n)]
        print("sweep", it, "residual 2^", m.bit_length() - PB, file=sys.stderr)
    print(nstr(mpf(X[0]) / S, 70))
else:  # Jacobi-PCG in fixed point; preconditioner z = r / p_i (integer division)
    b = [S] + [0]*(n-1)
    X = [0]*n; R = b[:]; Z = [R[i] // P[i] for i in range(n)]; D = Z[:]
    rz = sum(R[i]*Z[i] for i in range(n))
    for it in range(400):
        AD = matvec(D)
        dAd = sum(D[i]*AD[i] for i in range(n))
        if dAd == 0: break
        # alpha = rz / dAd as fixed point
        al = (rz << PB) // dAd
        X = [X[i] + (al*D[i] >> PB) for i in range(n)]
        R = [R[i] - (al*AD[i] >> PB) for i in range(n)]
        Z = [R[i] // P[i] for i in range(n)]
        rz2 = sum(R[i]*Z[i] for i in range(n))
        if max(abs(r) for r in R) < (1 << 20): break
        be = (rz2 << PB) // rz
        D = [Z[i] + (be*D[i] >> PB) for i in range(n)]
        rz = rz2
    print("cg its", it, file=sys.stderr)
    print(nstr(mpf(X[0]) / S, 70))
