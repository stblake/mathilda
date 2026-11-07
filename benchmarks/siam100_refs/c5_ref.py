"""c5 reference (method 1).  The benchmark's statement of Problem 5 is the
DISCRETE minimax problem on the M = 128 roots of unity z_k = exp(2 pi i k/M):

    E_M = min over complex c0..c3 of max_k |1/Gamma(z_k) - (c0 + c1 z_k + c2 z_k^2 + c3 z_k^3)|.

Method: (1) machine solve of the epigraph form (min t s.t. |e_k|^2 <= t^2) with
scipy SLSQP over the 8 real coefficients, (2) read off the active set,
(3) Newton (mpmath.findroot) at dps 80 on |e_k(c)|^2 = E^2 over the active
points, using the real-coefficient reduction (the sample set and 1/Gamma are
conjugation-symmetric, so the unique best approximation has real coefficients),
(4) certify: every inactive point has |e_k| < E, and the KKT multipliers of the
active set are positive.  Also prints the continuous answer's book value for
the gap."""
import numpy as np, mpmath as mp
from scipy.optimize import minimize
from scipy.special import rgamma
M = 128
th = 2*np.pi*np.arange(M)/M; z = np.exp(1j*th); fz = rgamma(z)
V = np.vander(z, 4, increasing=True)
def err(c): return fz - V @ (c[:4] + 1j*c[4:8])
x0 = np.zeros(9); x0[:4] = np.linalg.lstsq(V, fz, rcond=None)[0].real; x0[8] = 1
cons = {"type": "ineq", "fun": lambda x: x[8] - np.abs(err(x))}
r = minimize(lambda x: x[8], x0, constraints=[cons], bounds=[(None, None)]*8 + [(0, None)], method="SLSQP", options={"ftol": 1e-15, "maxiter": 1000})
E0 = r.x[8]; e = np.abs(err(r.x))
act = [k for k in range(M) if e[k] > E0 - 1e-7]
print("machine:", r.x, E0, "active:", act)
# reduce by symmetry: keep k in 0..M/2
actr = sorted(set(min(k, (M - k) % M) for k in act))
print("reduced active:", actr)
mp.mp.dps = 80
zs = [mp.expjpi(mp.mpf(2*k)/M) for k in range(M)]
fs = [mp.rgamma(w) for w in zs]
def e2(a, k):
    p = a[0] + a[1]*zs[k] + a[2]*zs[k]**2 + a[3]*zs[k]**3
    return abs(fs[k] - p)**2
assert len(actr) == 5, actr
def F(a0, a1, a2, a3, E):
    a = (a0, a1, a2, a3)
    return [e2(a, k) - E**2 for k in actr]
sol = mp.findroot(F, [mp.mpf(v) for v in r.x[:4]] + [mp.mpf(E0)])
a = sol[:4]; E = sol[4]
mx = max(mp.sqrt(e2(a, k)) for k in range(M) if min(k, M-k) not in actr)
print("max inactive / E:", mp.nstr(mx/E, 10))
# KKT: sum lam_k grad_a |e_k|^2 = 0, sum lam_k = 1, lam >= 0
J = mp.matrix(5, 5)
for row, k in enumerate(actr):
    for c in range(4):
        h = mp.mpf(10)**-30; ap = list(a); ap[c] += h; am = list(a); am[c] -= h
        J[c, row] = (e2(ap, k) - e2(am, k))/(2*h)
    J[4, row] = 1
lam = mp.lu_solve(J, mp.matrix([0, 0, 0, 0, 1]))
print("multipliers:", [mp.nstr(l, 6) for l in lam])
print("coeffs:", [mp.nstr(v, 20) for v in a])
print("E_128 =", mp.nstr(E, 70))
print("continuous (book App. B) = 0.2143352345904596394615264")
