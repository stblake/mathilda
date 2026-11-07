"""hc3 (holdout A, family C3) reference, method A: mpmath at dps 90.
sigma_max of each matrix by mpmath.eigsy on A^T A (Jacobi rotations), and the
closed forms for the two rank-one operators.  For the two infinite operators
with super-geometric/geometric decay (k = 3, 4) the section is taken at
n = 60 / 160, far past where it stops changing at 70 digits (also printed:
the benchmarked n = 40 / 100 sections).  Method B: hc3_ref.wl
(WL SingularValueList at precision 120)."""
import mpmath as mp
mp.mp.dps = 90
F = mp.mpf
fac = mp.factorial
E = {
    1: (lambda i, j: F(2)**-i * F(3)**-j, 100, 100),
    2: (lambda i, j: 1/(fac(i)*fac(j)), 30, 30),
    3: (lambda i, j: 1/fac(i + j - 1), 60, 60),
    4: (lambda i, j: 1/((i + j - 1)*F(2)**(i + j)), 160, 160),
    5: (lambda i, j: F(1)/(i + j - 1), 100, 100),
    6: (lambda i, j: F(min(i, j))/max(i, j), 100, 100),
    7: (lambda i, j: 1/F((i + j - 1)*(i + j - 2)//2 + j)**2, 100, 100),
    8: (lambda i, j: F(1)/(i + j*j), 120, 60),
    9: (lambda i, j: F(1)/(1 + abs(i - j))**2, 150, 150),
    10: (lambda i, j: mp.cos(i*j)/(i + j), 100, 100),
}


def smax(f, n, m):
    A = mp.matrix([[f(i, j) for j in range(1, m + 1)] for i in range(1, n + 1)])
    G = A.T * A
    ev = mp.eigsy(G, eigvals_only=True)
    return mp.sqrt(max(ev))


closed = {1: mp.sqrt(F(1)/24), 2: mp.besseli(0, 2) - 1}
for k in sorted(E):
    f, n, m = E[k]
    v = smax(f, n, m)
    extra = ""
    if k in closed:
        extra = " closed form agrees to %.1f digits" % float(-mp.log10(abs(v - closed[k])/closed[k]))
    if k == 3:
        extra = " (n=40 section differs by %s)" % mp.nstr(abs(smax(f, 40, 40) - v)/v, 3)
    if k == 4:
        extra = " (n=100 section differs by %s)" % mp.nstr(abs(smax(f, 100, 100) - v)/v, 3)
    print(k, mp.nstr(v, 70) + extra, flush=True)
