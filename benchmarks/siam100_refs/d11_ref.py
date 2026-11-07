"""Reference for d11 (six masses), method 2 (independent of d11_ref.wl's
NDSolve): mpmath.odefun (Taylor-series integrator) at dps 75, shot in
(theta, T) by mpmath.findroot (multidimensional Newton), started on the
shortest branch found by d11_scan.py."""
import time
from mpmath import mp, mpf, cos, sin, odefun, findroot, nstr
mp.dps = 75
M = [(2, -1), (2, 0), (2, 1), (3, -1), (3, 0), (3, 1)]
def rhs(t, u):
    x, y, vx, vy = u
    ax = ay = mpf(0)
    for (a, b) in M:
        dx, dy = a - x, b - y
        r3 = (dx * dx + dy * dy) ** mpf(1.5)
        ax += dx / r3; ay += dy / r3
    return [vx, vy, ax, ay]
cache = {}
def traj(th):
    if th not in cache:
        cache[th] = odefun(rhs, 0, [mpf(0), mpf(0), cos(th), sin(th)])
    return cache[th]
def F(th, T):
    u = traj(th)(T)
    return [u[0] - 4, u[1] - 1]
t0 = time.time()
sol = findroot(F, (mpf("-0.274881254"), mpf("1.9160969183")), tol=mpf(10) ** -140)
print("theta", nstr(sol[0], 70))
print("T    ", nstr(sol[1], 70))
print("residual", [nstr(v, 5) for v in F(sol[0], sol[1])], "secs", round(time.time() - t0))
