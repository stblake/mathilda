"""c3 reference for the reduced rows (method 2): sigma_max of the n x n
section A_n, a_ij = 1/((i+j-1)(i+j-2)/2 + j), in mpmath:
  n = 100: mpmath.svd_r singular values at dps 80;
  n = 500: power iteration on A^T A with exact-rational entries rounded to
           dps 90 (sigma2/sigma1 ~ 0.11), until the Rayleigh quotient is stable.
Compare with c3_trunc.m (WL, 160-digit power iteration: method 1)."""
import mpmath as mp
mp.mp.dps = 80
def A(n): return [[mp.mpf(1)/((i+j-1)*(i+j-2)//2 + j) for j in range(1, n+1)] for i in range(1, n+1)]
s = mp.svd_r(mp.matrix(A(100)), compute_uv=False)
print("n=100", mp.nstr(max(s), 70))
mp.mp.dps = 90
n = 500; a = A(n); at = [list(c) for c in zip(*a)]
mv = lambda m, v: [mp.fsum(x*y for x, y in zip(row, v)) for row in m]
v = [mp.mpf(1)]*n; lam0 = 0
for k in range(200):
    w = mv(at, mv(a, v)); lam = mp.sqrt(mp.fsum(x*x for x in w)/mp.fsum(x*x for x in v))
    v = w
    if abs(lam - lam0) < mp.mpf(10)**-85: break
    lam0 = lam
print("n=500", mp.nstr(mp.sqrt(lam), 70), "iters", k)
