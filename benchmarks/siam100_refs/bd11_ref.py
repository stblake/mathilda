"""Reference for holdout-B family d11 (experiment 246): shooting in two (or three)
unknowns with an ODE solve inside a vector root finder.

For every problem two independent methods:
  A  closed form / analytic reduction where one exists, else mpmath.odefun
     (Taylor) inside mpmath.findroot;
  B  mpmath.odefun shot (when A is analytic) or the Gragg-Bulirsch-Stoer
     extrapolation integrator of bhold_ode_util.py shot by findroot.
Usage: bd11_ref.py k [dps]   (prints A, B and their agreement for problem k)
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, cos, sin, exp, sqrt, asin, acosh, ellipk, odefun, findroot, nstr, pi, atan2, mpc
from bhold_ode_util import bs_solve

k = int(sys.argv[1]); mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 75
TOL = mpf(10) ** (-(mp.dps - 10))


def shoot(rhs, init, cond, guess, integ):
    """findroot over z of cond(z, y(T)) where y solves rhs from init(z)."""
    def F(*z):
        y0, T = init(z)
        if integ == "taylor":
            yT = odefun(rhs, 0, y0)(T)
        else:
            yT = bs_solve(rhs, 0, y0, T)
        return cond(z, yT)
    return findroot(F, [mpf(g) for g in guess], tol=TOL ** 2 if False else mpf(10) ** (-2 * mp.dps + 20))


def planar(rhs_acc, speed, start, target, guess, integ):
    rhs = lambda t, u: [u[2], u[3]] + rhs_acc(u[0], u[1], u[2], u[3])
    init = lambda z: ([mpf(start[0]), mpf(start[1]), speed * cos(z[0]), speed * sin(z[0])], z[1])
    cond = lambda z, y: [y[0] - target[0], y[1] - target[1]]
    return shoot(rhs, init, cond, guess, integ)


out = {}
t0 = time.time()
if k == 1:   # linear drag projectile, v = 2, k = 1/2, g = 1, target (1, 1/2)
    def cf(th, T):
        E = 1 - exp(-T / 2)
        return [4 * cos(th) * E - 1, (4 * sin(th) + 4) * E - 2 * T - mpf(1) / 2]
    A = findroot(cf, [mpf("0.45"), mpf("0.6")])
    B = planar(lambda x, y, vx, vy: [-vx / 2, -1 - vy / 2], 2, (0, 0), (1, mpf(1) / 2), ("0.45", "0.6"), "taylor")
    out = {"A closed form x(t),y(t) + findroot": A, "B odefun shot": B}
elif k == 2:  # isotropic harmonic, start (1,0), speed 1, target (0, 1/2)
    A = asin(sqrt(10) / 4)
    B = planar(lambda x, y, vx, vy: [-x, -y], 1, (1, 0), (0, mpf(1) / 2), ("2.47", "0.91"), "taylor")
    out = {"A T = asin(sqrt(10)/4)": A, "B odefun shot": B}
elif k == 3:  # Kepler, start (1,0), speed 11/10, target (-1/2, 1)
    v = mpf(11) / 10; a = 1 / (2 - v * v); r1 = mpf(1); r2 = sqrt(mpf(5)) / 2
    c = sqrt((mpf(-1) / 2 - 1) ** 2 + 1); s = (r1 + r2 + c) / 2
    al = 2 * asin(sqrt(s / (2 * a))); be = 2 * asin(sqrt((s - c) / (2 * a)))
    cands = [a ** 1.5 * ((al - sin(al)) - (be - sin(be))), a ** 1.5 * ((2 * pi - al + sin(al)) - (be - sin(be)))]
    A = min(cands)
    B = planar(lambda x, y, vx, vy: [-x / (x * x + y * y) ** mpf(1.5), -y / (x * x + y * y) ** mpf(1.5)],
               v, (1, 0), (mpf(-1) / 2, 1), ("1.71", "1.84"), "taylor")
    out = {"A Lambert/Lagrange time equation (a = 100/79), candidates " + str([nstr(x, 12) for x in cands]): A,
           "B odefun shot": B}
elif k == 4:  # Bratu y'' = -lam e^y, y(0)=0, y'(0)=s; y(1)=0, y(1/2)=1/2  -> lam
    A = 8 * acosh(exp(mpf(1) / 4)) ** 2 * exp(mpf(-1) / 2)
    rhs = lambda t, u: [u[1], -u[2] * exp(u[0]), mpf(0)]
    def F(lam, s):
        f = odefun(rhs, 0, [mpf(0), s, lam])
        return [f(1)[0], f(mpf(1) / 2)[0] - mpf(1) / 2]
    B = findroot(F, [mpf("2.7"), mpf("1.9")])
    out = {"A 8 acosh(e^(1/4))^2 e^(-1/2) (Liouville/Gelfand closed form)": A, "B odefun shot": B}
elif k == 5:  # clamped beam y'''' + 50 (1+x) y = 1 -> y''(0)
    from fractions import Fraction as Q
    # exact power series: y = sum c_n x^n, c0=c1=0, c2=p/2, c3=q/6; linear in (1, p, q)
    N = 400
    def series(c2, c3, rhs1):
        c = [Q(0)] * (N + 5); c[2] = c2; c[3] = c3
        for n in range(0, N):
            # (n+4)(n+3)(n+2)(n+1) c_{n+4} = [n==0]*rhs1 - 50 (c_n + c_{n-1})
            r = (Q(rhs1) if n == 0 else Q(0)) - 50 * (c[n] + (c[n - 1] if n >= 1 else 0))
            c[n + 4] = r / ((n + 4) * (n + 3) * (n + 2) * (n + 1))
        return c
    def ends(c):
        return sum(c), sum(n * cn for n, cn in enumerate(c))
    y0 = ends(series(Q(0), Q(0), 1)); yp = ends(series(Q(1, 2), Q(0), 0)); yq = ends(series(Q(0), Q(1, 6), 0))
    # solve y0 + p yp + q yq = 0 exactly
    det = yp[0] * yq[1] - yq[0] * yp[1]
    p = (-y0[0] * yq[1] + yq[0] * y0[1]) / det
    A = mpf(p.numerator) / p.denominator
    rhs = lambda t, u: [u[1], u[2], u[3], 1 - 50 * (1 + t) * u[0]]
    init = lambda z: ([mpf(0), mpf(0), z[0], z[1]], 1)
    B = shoot(rhs, init, lambda z, y: [y[0], y[1]], ("0.08", "-0.5"), "taylor")
    out = {"A exact rational power series (400 terms) + 2x2 solve": A, "B odefun shot": B}
elif k == 6:  # pendulum th''=-sin th, th(0)=0, th'(0)=w; th(T)=2, th'(T)=0 -> T
    A = ellipk(sin(1) ** 2)
    rhs = lambda t, u: [u[1], -sin(u[0])]
    init = lambda z: ([mpf(0), z[0]], z[1])
    B = shoot(rhs, init, lambda z, y: [y[0] - 2, y[1]], ("1.7", "2.1"), "taylor")
    out = {"A K(m = sin(1)^2)": A, "B odefun shot": B}
elif k == 7:  # Henon-Heiles, speed 1/2, target (3/10, 1/5)
    acc = lambda x, y, vx, vy: [-x - 2 * x * y, -y - x * x + y * y]
    A = planar(acc, mpf(1) / 2, (0, 0), (mpf(3) / 10, mpf(1) / 5), ("0.5836", "0.8305"), "taylor")
    B = planar(acc, mpf(1) / 2, (0, 0), (mpf(3) / 10, mpf(1) / 5), ("0.5836", "0.8305"), "bs")
    out = {"A odefun shot": A, "B Bulirsch-Stoer shot": B}
elif k == 8:  # two unit masses at (+-1, 0), start (0,-2), speed 4/5, target (1/2, 3/2)
    def acc(x, y, vx, vy):
        ax = ay = mpf(0)
        for (p, q) in ((-1, 0), (1, 0)):
            dx, dy = p - x, q - y; r2 = dx * dx + dy * dy; r3 = r2 * sqrt(r2)
            ax += dx / r3; ay += dy / r3
        return [ax, ay]
    A = planar(acc, mpf(4) / 5, (0, -2), (mpf(1) / 2, mpf(3) / 2), ("1.50", "2.75"), "taylor")
    B = planar(acc, mpf(4) / 5, (0, -2), (mpf(1) / 2, mpf(3) / 2), ("1.50", "2.75"), "bs")
    out = {"A odefun shot": A, "B Bulirsch-Stoer shot": B}
elif k == 9:  # 3-D x''=-x, y''=-4y, z''=-9z, from 0 speed 1, target (3/10, 1/5, 1/10)
    g = lambda t: (mpf(3) / 10 / sin(t)) ** 2 + (mpf(2) / 5 / sin(2 * t)) ** 2 + (mpf(3) / 10 / sin(3 * t)) ** 2 - 1
    A = findroot(g, mpf("0.37"))
    a = mpf(3) / 10 / sin(A); b = mpf(2) / 5 / sin(2 * A); c = mpf(3) / 10 / sin(3 * A)
    th0, ph0 = atan2(b, a), asin(c)
    rhs = lambda t, u: [u[3], u[4], u[5], -u[0], -4 * u[1], -9 * u[2]]
    init = lambda z: ([mpf(0)] * 3 + [cos(z[1]) * cos(z[0]), cos(z[1]) * sin(z[0]), sin(z[1])], z[2])
    cond = lambda z, y: [y[0] - mpf(3) / 10, y[1] - mpf(1) / 5, y[2] - mpf(1) / 10]
    B = shoot(rhs, init, cond, ("0.6", "0.25", "0.38"), "taylor")
    print("start", nstr(th0, 10), nstr(ph0, 10))
    out = {"A 1-D trig equation (a/sin t)^2+... = 1, findroot": A, "B odefun shot (theta,phi,T)": B}
elif k == 10:  # Arenstorf orbit half period: x0 = 0.994, y'(0)=v; y(T)=0, x'(T)=0
    mu = mpf(12277471) / 10 ** 9; m = 1 - mu
    def rhs(t, u):
        x, y, vx, vy = u
        d1 = (x + mu) ** 2 + y * y; d1 = d1 * sqrt(d1); d2 = (x - m) ** 2 + y * y; d2 = d2 * sqrt(d2)
        return [vx, vy, x + 2 * vy - m * (x + mu) / d1 - mu * (x - m) / d2, y - 2 * vx - m * y / d1 - mu * y / d2]
    init = lambda z: ([mpf("0.994"), mpf(0), mpf(0), z[0]], z[1])
    cond = lambda z, y: [y[1], y[2]]
    which = sys.argv[3] if len(sys.argv) > 3 else "taylor"
    A = shoot(rhs, init, cond, ("-2.0015851063790825", "8.532608280078981"), which)
    out = {which + " shot (published full period 17.0652165601579625588917206249, Hairer-Norsett-Wanner)": A}
for kk, vv in out.items():
    print(kk); print("   ", nstr(vv, 70) if not isinstance(vv, mp.matrix) else [nstr(x, 70) for x in vv])
print("secs", round(time.time() - t0))
