"""Reference A for SIAM100 App. D problem 19 (exp 194): the finite proxy.

Problem 19 (book): -div(c grad u) = 1 on the unit square, u = 0 on the
boundary; inf of Int u over measurable c in {1, 100}, each on half the area.
No value is known.  Relaxation: for ONE load in a scalar problem, the optimal
microstructure is a rank-1 laminate aligned with the flux, whose conductivity
along the flux is the ARITHMETIC mean -- so the relaxed problem is the convex
"variable thickness" problem  min { J(c) : 1 <= c <= 100, mean(c) = 101/2 },
J(c) = Int u = compliance, convex in c.

THE PROXY (what exp 194 asks): n x n uniform squares (n = 8), c constant on
each, bilinear (Q1) finite elements for u, J_n(c) = F . U with K(c) U = F,
element stiffness c * (1/6)[[4,-1,-2,-1],[-1,4,-1,-2],[-2,-1,4,-1],[-1,-2,-1,4]]
(corners counter-clockwise), F_i = h^2 at every interior node.  Minimise J_n
over c in [1,100]^(n^2) with sum(c) = 50.5 n^2.  Convex => unique value.

Method A (this file): scipy SLSQP locates the optimum and its active set
(c_e = 1, c_e = 100, or interior); then the KKT system on the interior set
    U^T K_e U = mu  (e interior),  sum c = 50.5 n^2,  K(c) U = F
is solved by mpmath.findroot (Newton) at 80 digits -- reduced to one unknown
per orbit of the square's symmetry group D4 (the problem is D4-invariant and
convex, so the group average of any optimum is an optimum), and the KKT sign
conditions are verified at that precision (U^T K_e U >= mu where c_e = 100,
<= mu where c_e = 1).
Usage: python3.11 d19_ref.py [n]
"""
import sys
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
h = 1.0 / n
Kq = np.array([[4, -1, -2, -1], [-1, 4, -1, -2], [-2, -1, 4, -1], [-1, -2, -1, 4]]) / 6.0


def node(i, j):           # interior node index or -1
    if 1 <= i <= n - 1 and 1 <= j <= n - 1:
        return (i - 1) * (n - 1) + (j - 1)
    return -1


elems = [[node(a, b), node(a + 1, b), node(a + 1, b + 1), node(a, b + 1)]
         for a in range(n) for b in range(n)]
N = (n - 1) ** 2


def assemble(c, zero, one, six):
    K = [[zero] * N for _ in range(N)]
    for e, idx in enumerate(elems):
        for r in range(4):
            if idx[r] < 0:
                continue
            for s in range(4):
                if idx[s] < 0:
                    continue
                K[idx[r]][idx[s]] += c[e] * KQ[r][s] / six
    return K


KQ = [[4, -1, -2, -1], [-1, 4, -1, -2], [-2, -1, 4, -1], [-1, -2, -1, 4]]


def J_np(c):
    K = np.array(assemble(c, 0.0, 1.0, 6.0))
    F = np.full(N, h * h)
    U = np.linalg.solve(K, F)
    return F @ U, U


def grad_np(c):
    _, U = J_np(c)
    g = np.zeros(len(c))
    for e, idx in enumerate(elems):
        ue = np.array([U[k] if k >= 0 else 0.0 for k in idx])
        g[e] = -ue @ Kq @ ue
    return g


c0 = np.full(n * n, 50.5)
res = minimize(lambda c: J_np(c)[0], c0, jac=grad_np, method="SLSQP",
               bounds=[(1, 100)] * (n * n),
               constraints=[{"type": "eq", "fun": lambda c: c.sum() - 50.5 * n * n,
                             "jac": lambda c: np.ones_like(c)}],
               options={"ftol": 1e-16, "maxiter": 2000})
print("SLSQP:", res.message, repr(res.fun), flush=True)
c = res.x
lo = [e for e in range(n * n) if c[e] < 1 + 1e-6]
hi = [e for e in range(n * n) if c[e] > 100 - 1e-6]
free = [e for e in range(n * n) if e not in lo and e not in hi]
print("active: c=1:", len(lo), " c=100:", len(hi), " free:", len(free))

mp.mp.dps = 80
F = mp.matrix([mp.mpf(1) / n ** 2] * N)


def solveU(cv):
    K = mp.matrix(assemble(cv, mp.mpf(0), 1, mp.mpf(6)))
    return mp.lu_solve(K, F)


def energy(U, e):
    ue = [U[k] if k >= 0 else mp.mpf(0) for k in elems[e]]
    return sum(ue[r] * KQ[r][s] * ue[s] for r in range(4) for s in range(4)) / 6


def orbit(e):
    a, b = divmod(e, n)
    imgs = set()
    for (p, q) in [(a, b), (b, a)]:
        for pp in (p, n - 1 - p):
            for qq in (q, n - 1 - q):
                imgs.add(pp * n + qq)
    return min(imgs)


orbits = sorted(set(orbit(e) for e in free))
members = {o: [e for e in free if orbit(e) == o] for o in orbits}
print("free orbits:", len(orbits), [len(members[o]) for o in orbits], flush=True)


def full_c(xorb):
    cv = [mp.mpf(1)] * (n * n)
    for e in hi:
        cv[e] = mp.mpf(100)
    for k, o in enumerate(orbits):
        for e in members[o]:
            cv[e] = xorb[k]
    return cv


def kkt(*args):
    xs, mu = list(args[:-1]), args[-1]
    U = solveU(full_c(xs))
    # scaled to O(1): element energies are ~1e-7, mu is carried as mu * SCALE
    eqs = [energy(U, o) * SCALE - mu for o in orbits]
    eqs.append((sum(full_c(xs)) - mp.mpf(101) / 2 * n * n) / (n * n))
    return eqs


g = grad_np(c)
mu0 = float(np.mean([-g[e] for e in free])) if free else 0.0
SCALE = mp.mpf(10) ** 7
x0 = [mp.mpf(float(np.mean([c[e] for e in members[o]]))) for o in orbits] + [mp.mpf(mu0) * SCALE]
sol = mp.findroot(kkt, x0, tol=mp.mpf(10) ** -140, maxsteps=50)
sol = [sol[i] for i in range(len(x0))]   # (mp.matrix does not support negative indices)
xs, mu = sol[:-1], sol[-1] / SCALE
print("KKT residual:", mp.nstr(max(abs(v) for v in kkt(*sol)), 5))
cv = full_c(xs)
U = solveU(cv)
Jv = mp.fdot(F, U)
ok_hi = all(energy(U, e) >= mu for e in hi)
ok_lo = all(energy(U, e) <= mu for e in lo)
ok_free = all(1 < x < 100 for x in xs) and all(abs(energy(U, e) - mu) < mp.mpf(10) ** -70 for e in free)
print("KKT signs: hi", ok_hi, " lo", ok_lo, " free in (1,100)", ok_free, " mu", mp.nstr(mu, 20))
print("J_n =", mp.nstr(Jv, 70))
# hand the optimal design to the independent Wolfram-side certificate
# (d19_ref.wl): one value per cell, element order e = a*n + b.
with open("d19_cstar_n%d.txt" % n, "w") as fh:
    fh.write("\n".join(mp.nstr(v, 78) for v in cv) + "\n")
