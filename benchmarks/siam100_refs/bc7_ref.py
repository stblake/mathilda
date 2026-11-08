"""Holdout B, family C7 (experiment 232): references for bc7.1 .. bc7.10.
Matrices: bc7_mats.py (built exactly as in the benchmark halves).

Method A (all ten, `refine`): mixed-precision iterative refinement with the
residual b - A x computed EXACTLY in Python integers (x held as an integer
scaled by 2^PB); the correction solved in doubles by scipy splu (or
Jacobi-PCG for n = 20000).  Each sweep gains ~12-15 digits; run to a zero
residual at 2^-PB.

Method B (independent, per problem):
  bc7.1   the exact finite continued fraction 1/(p1 - 1/(p2 - ...)) in Fractions
  bc7.2/8/9  the eigen-decomposition (discrete sine transform) closed-form sums
  bc7.3/6 Jacobi-PCG in exact fixed-point integer arithmetic (no doubles)
  bc7.4/5 PARI/gp matsolve (dense Gaussian elimination, realprecision 80)
  bc7.7   the character sum over Z_p^*  (1/(p-1)) sum_k 1/(7 + sum_{h=2,3,5} 2cos(2 pi k e_h/(p-1))),
          h = g^{e_h} for a primitive root g (the matrix is a convolution on the cyclic group)
  bc7.10  the exact level recursion (Schur complement up the tree) in Fractions
Usage: bc7_ref.py k   (prints A, B and agreement)
"""
import sys, subprocess
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from fractions import Fraction
import mpmath
from mpmath import mp, mpf, nstr
from bc7_mats import build, primes

mp.dps = 80
PB = 330
S = 1 << PB
k = int(sys.argv[1])
n, E, b, idx = build(k)
rows = [[] for _ in range(n)]
for (i, j), v in E.items():
    rows[i].append((j, v))


def matvec(x):
    return [sum(v * x[j] for j, v in rows[i]) for i in range(n)]


def refine():
    ii, jj, vv = zip(*[(i, j, float(v)) for (i, j), v in E.items()])
    A = sp.csc_matrix((vv, (ii, jj)), shape=(n, n))
    if n > 10000 and k in (3,):
        M = sp.diags([1.0 / A[i, i] for i in range(n)])
        solve = lambda r: sla.cg(A, r, rtol=1e-15, atol=0.0, M=M, maxiter=5000)[0]
    else:
        lu = sla.splu(A); solve = lu.solve
    X = [0] * n
    for it in range(12):
        AX = matvec(X)
        R = [b[i] * S - AX[i] for i in range(n)]
        m = max(abs(r) for r in R)
        print("  sweep", it, "residual 2^", m.bit_length() - PB if m else None, file=sys.stderr)
        if m == 0 or m.bit_length() < 8: break
        e = m.bit_length() - 60
        rf = np.array([float(Fraction(r, 1 << e)) if e > 0 else float(r) for r in R])
        d = solve(rf)
        X = [X[i] + int(round(d[i] * 2.0 ** 50)) * (1 << e) // (1 << 50) if e >= 0
             else X[i] + int(round(d[i] * 2.0 ** e)) for i in range(n)]
    return mpf(X[idx]) / S


def fixed_pcg():
    D0 = [dict(rows[i])[i] for i in range(n)]
    X = [0] * n; R = [bi * S for bi in b]
    Z = [R[i] // D0[i] for i in range(n)]; Dv = Z[:]
    rz = sum(R[i] * Z[i] for i in range(n))
    for it in range(3000):
        AD = matvec(Dv)
        dAd = sum(Dv[i] * AD[i] for i in range(n))
        al = (rz << PB) // dAd
        X = [X[i] + (al * Dv[i] >> PB) for i in range(n)]
        R = [R[i] - (al * AD[i] >> PB) for i in range(n)]
        if max(abs(r) for r in R) < (1 << 24): break
        Z = [R[i] // D0[i] for i in range(n)]
        rz2 = sum(R[i] * Z[i] for i in range(n))
        be = (rz2 << PB) // rz
        Dv = [Z[i] + (be * Dv[i] >> PB) for i in range(n)]
        rz = rz2
    print("  pcg its", it, file=sys.stderr)
    return mpf(X[idx]) / S


def sine_sum(N, d, node, shift):
    """(A^-1)_{node,node}, A = sum of d copies of T_N (+shift I); node = tuple of 1-based coords."""
    h = mpf(2) / (N + 1)
    lam = [4 * mpmath.sin(j * mpmath.pi / (2 * (N + 1))) ** 2 for j in range(1, N + 1)]
    w = [[h * mpmath.sin(j * mpmath.pi * c / (N + 1)) ** 2 for j in range(1, N + 1)] for c in node]
    if d == 2:
        return mpmath.fsum(w[0][a] * w[1][bb] / (lam[a] + lam[bb] + shift)
                           for a in range(N) for bb in range(N))
    return mpmath.fsum(w[0][a] * w[1][bb] * w[2][c] / (lam[a] + lam[bb] + lam[c] + shift)
                       for a in range(N) for bb in range(N) for c in range(N))


def gp_dense():
    P = primes(n)
    if k == 4:
        offs = [j * j for j in range(1, 45) if j * j < n]; sub, sup = 1, 1
    else:
        offs = [2 ** j for j in range(11)]; sub, sup = -1, 1
    prog = ("default(realprecision, 80); n=%d; M=matrix(n,n); for(i=1,n,M[i,i]=prime(i)*1.);"
            "foreach(%s, o, for(i=1,n-o, M[i,i+o]=%d; M[i+o,i]=%d));"
            "v=vectorv(n); v[1]=1; x=matsolve(M,v); printf(\"%%.75g\\n\", x[1]); quit\n") % (n, str(offs), sup, sub)
    out = subprocess.run(["/usr/local/bin/gp", "-q", "-s", "400000000", "--default", "parisizemax=6000000000"], input=prog, capture_output=True, text=True).stdout
    return mpf(out.strip().split()[-1])


a = refine()
if k == 1:
    P = primes(n); g = Fraction(P[-1])
    for p in reversed(P[:-1]):
        g = p - 1 / g
    bval = mpmath.mpf(g.denominator) / g.numerator
elif k == 2:
    bval = sine_sum(100, 2, (50, 50), 1)
elif k in (3, 6):
    bval = fixed_pcg()
elif k in (4, 5):
    bval = gp_dense()
elif k == 7:
    p = 2003
    from sympy import primitive_root, discrete_log
    g = primitive_root(p); es = [discrete_log(p, h, g) for h in (2, 3, 5)]
    tp = 2 * mpmath.pi / (p - 1)
    bval = mpmath.fsum(1 / (7 + sum(2 * mpmath.cos(tp * ((j * e) % (p - 1))) for e in es)) for j in range(p - 1)) / (p - 1)
elif k == 8:
    bval = sine_sum(30, 3, (1, 1, 1), 0)
elif k == 9:
    bval = sine_sum(150, 2, (75, 75), mpf(-1) / 2)
elif k == 10:
    g = Fraction(2)                      # leaves: degree 1, +1
    for level in range(9):               # levels 9 .. 1 (internal, degree 3, +1)
        g = 4 - 2 / g
    g = 3 - 2 / g                         # root: degree 2, +1
    bval = mpmath.mpf(g.denominator) / g.numerator
print("bc7.%d A" % k, nstr(a, 72))
print("bc7.%d B" % k, nstr(bval, 72))
print("agree", nstr(-mpmath.log10(abs(a - bval) / abs(a)), 3))
