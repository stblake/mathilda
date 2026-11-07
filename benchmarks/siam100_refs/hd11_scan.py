"""Holdout A, family d11 (experiment 216): branch scan for the shooting
problems without a closed form (hd11.7, hd11.8, hd11.10).  For each launch
angle on a grid, integrate (DOP853, rtol 1e-10) to t = Tmax with dense output,
locate the local minima in t of the distance to the target, refine every near
hit (distance < 0.05) with fsolve in (theta, T), and print the distinct
solutions by increasing T.  The SHORTEST is the problem's answer."""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# k: (masses, target[, start, speed])
CASES = {5: ([(0, 0)], (-0.5, 0.5), (1, 0), 1.0),
         6: ([(0, 0)], (0, 1.5), (1, 0), 1.2),
         7: ([(2, -1), (2, 1)], (4, 0)),
         8: ([(1, 1), (2, -1), (3, 1)], (4, 0)),
         10: ([(2, -1), (2, 0), (2, 1), (3, -1), (3, 0), (3, 1)], (1, 2))}

def make(ms):
    M = np.array(ms, float)
    def rhs(t, u):
        d = M - u[:2]
        a = (d / ((d * d).sum(1) ** 1.5)[:, None]).sum(0)
        return [u[2], u[3], a[0], a[1]]
    return rhs

START, SPEED = (0, 0), 1.0
def traj(rhs, th, T, dense=False):
    return solve_ivp(rhs, [0, T], [START[0], START[1], SPEED * np.cos(th), SPEED * np.sin(th)], method="DOP853",
                     rtol=1e-10, atol=1e-10, dense_output=dense)

def scan(k, n=1440, Tmax=4.0):
    global START, SPEED
    ms, tgt = CASES[k][:2]
    START, SPEED = (CASES[k][2], CASES[k][3]) if len(CASES[k]) > 2 else ((0, 0), 1.0)
    rhs = make(ms)
    sols = []
    tt = np.linspace(0, Tmax, 4001)
    for th in np.linspace(-np.pi, np.pi, n, endpoint=False):
        s = traj(rhs, th, Tmax, True)
        te = tt[tt <= s.t[-1]]
        if len(te) < 3: continue
        p = s.sol(te)
        d = np.hypot(p[0] - tgt[0], p[1] - tgt[1])
        for i in range(1, len(d) - 1):
            if d[i] < d[i - 1] and d[i] <= d[i + 1] and d[i] < 0.05:
                def F(z):
                    e = traj(rhs, z[0], z[1])
                    return [e.y[0, -1] - tgt[0], e.y[1, -1] - tgt[1]]
                z, info, ier, _ = fsolve(F, [th, te[i]], full_output=True, xtol=1e-12)
                if ier == 1 and z[1] > 0 and max(abs(np.array(F(z)))) < 1e-8:
                    zz = (np.mod(z[0] + np.pi, 2 * np.pi) - np.pi, z[1])
                    if not any(abs(zz[1] - q[1]) < 1e-7 for q in sols):
                        sols.append(zz)
    for th, T in sorted(sols, key=lambda q: q[1]):
        print("hd11.%d  theta % .12f  T %.12f" % (k, th, T))

for k in map(int, sys.argv[1:] or ["7", "8", "10"]):
    scan(k)
