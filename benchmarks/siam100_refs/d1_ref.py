"""D1 reference: sum of 1/n over n >= 1 whose decimal digits avoid "42".
Schmelzer & Baillie (Amer. Math. Monthly 115 (2008) 525-540) power-sum recursion.
States: s0 = last digit not 4, s1 = last digit 4.  psi_i(k, j) = sum of n^-j over
the k-digit admissible n in state i.  Appending a digit d to n:
  (10 n + d)^-j = sum_w (-1)^w C(j+w-1, w) d^w 10^(-j-w) n^(-j-w),
so psi(k+1, .) = M psi(k, .) with a CONSTANT matrix M (moments truncated at J).
Explicit enumeration for k <= K0, then the whole tail in closed form:
  S = sum_{k<=K0} psi(k,1) + e^T (I - M)^-1 M psi(K0, .).
Usage: d1_ref.py K0 J dps"""
import sys
from mpmath import mp, mpf, matrix, binomial, lu_solve, nstr
K0, J, dps = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
mp.dps = dps
psi = [[[mpf(0)]*(J+1) for _ in range(2)] for _ in range(K0+1)]
S = mpf(0)
for n in range(1, 10**K0):
    s = str(n)
    if "42" in s: continue
    k = len(s); st = 1 if s[-1] == "4" else 0
    x = mpf(1)/n; p = x
    for j in range(1, J+1):
        psi[k][st][j] += p; p *= x
    S += x
# transition list (from state, digit, to state)
T = []
for d in range(10):
    T.append((0, d, 1 if d == 4 else 0))
    if d != 2: T.append((1, d, 1 if d == 4 else 0))
N = 2*J
M = matrix(N, N)
idx = lambda st, j: st*J + (j-1)
for (i, d, l) in T:
    for j in range(1, J+1):
        for w in range(0, J-j+1):
            c = (-1)**w * binomial(j+w-1, w) * mpf(d)**w / mpf(10)**(j+w)
            M[idx(l, j), idx(i, j+w)] += c
v = matrix([psi[K0][st][j] for st in range(2) for j in range(1, J+1)])
Mv = M*v
I = matrix(N, N)
for i in range(N): I[i, i] = 1
y = lu_solve(I - M, Mv)
S += y[idx(0, 1)] + y[idx(1, 1)]
print(nstr(S, dps-10))
