"""hc5 (holdout A, family C5) reference, method A.  Every problem is a
DISCRETE minimax (best uniform) approximation by a polynomial of degree n:

    E = min over coefficients c of max_k |f(x_k) - sum_j c_j x_k^j|

Real problems (k = 1..6, real coefficients, real sample points): the LP
  min E s.t. -E <= f_k - p(x_k) <= E  (scipy linprog / HiGHS) gives the active
  set; on it the alternation system f_k - p(x_k) = s_k E (n+2 equations,
  linear in c, E) is solved by mpmath.lu_solve at dps 90; certified: every
  inactive |e_k| < E, and the LP dual multipliers are positive.
Complex problems (k = 7..10, complex coefficients, points on the unit circle):
  SLSQP on the epigraph form gives the active set; the data are
  conjugation-symmetric, so the unique optimum has real coefficients and the
  KKT system |e_k|^2 = E^2 (active k modulo conjugation),
  sum_k lam_k grad_a|e_k|^2 = 0, sum lam_k = 1 in the n+1 real coefficients,
  E and the lam_k is solved by mpmath.findroot at dps 90; certified: lam_k > 0 and every inactive |e_k| < E (convex problem =>
  global).  Closed form: k = 1 (E = 1/8: the continuous extremal points 0,
  1/4, 1 are samples).
Method B: hc5_ref.wl (exact-rational simplex LinearProgramming on the data
rounded to 95 digits for the real problems; FindRoot on the KKT system at
WorkingPrecision 120 for the complex ones)."""
import numpy as np
import mpmath as mp
from scipy.optimize import linprog, minimize
from scipy.special import rgamma

mp.mp.dps = 90
F = mp.mpf


def real_probs():
    P = {}
    P[1] = (lambda x: mp.sqrt(x), [F(k)/128 for k in range(129)], 1)
    P[2] = (lambda x: mp.exp(x), [mp.cos(k*mp.pi/64) for k in range(65)], 3)
    P[3] = (lambda x: abs(x), [F(-1) + F(k)/64 for k in range(129)], 4)
    P[4] = (lambda x: mp.rgamma(x), [F(2*k)/100 for k in range(101)], 3)
    P[5] = (lambda x: mp.atan(5*x), [mp.cos(k*mp.pi/100) for k in range(101)], 5)
    P[6] = (lambda x: mp.log(1 + x), [F(k)/64 for k in range(65)], 2)
    return P


def complex_probs():
    P = {}
    z = lambda M: [mp.expjpi(F(2*k)/M) for k in range(M)]
    P[7] = (lambda w: mp.rgamma(w), z(128), 2)
    P[8] = (lambda w: mp.gamma(w + 2), z(128), 2)
    P[9] = (lambda w: mp.log(2 + w), z(64), 2)
    P[10] = (lambda w: mp.rgamma(w), z(128), 4)
    return P


def solve_real(f, xs, n):
    fx = [f(x) for x in xs]
    X = np.array([float(x) for x in xs]); Fv = np.array([float(v) for v in fx])
    V = np.vander(X, n + 1, increasing=True)
    m = len(X)
    # vars c_0..c_n, E ; minimise E ; f - Vc <= E ; -(f - Vc) <= E
    A = np.vstack([np.hstack([-V, -np.ones((m, 1))]), np.hstack([V, -np.ones((m, 1))])])
    b = np.concatenate([-Fv, Fv])
    r = linprog(np.r_[np.zeros(n + 1), 1], A_ub=A, b_ub=b, bounds=[(None, None)]*(n + 2), method="highs")
    c, E0 = r.x[:n + 1], r.x[n + 1]
    e = Fv - V @ c
    act = [k for k in range(m) if abs(abs(e[k]) - E0) < 1e-9*max(E0, 1e-300) + 1e-13]
    assert len(act) >= n + 2, (len(act), act)
    # exact alternation system: sum_j c_j x_k^j + s_k E = f_k (square when
    # n + 2 points are active; an even/odd f on a symmetric set has n + 3
    # active points and a consistent overdetermined system: least squares,
    # residual checked)
    Amp = mp.matrix(len(act), n + 2); rhs = mp.matrix(len(act), 1)
    for row, k in enumerate(act):
        for j in range(n + 1):
            Amp[row, j] = xs[k]**j
        Amp[row, n + 1] = 1 if e[k] > 0 else -1
        rhs[row] = fx[k]
    sol, res = mp.qr_solve(Amp, rhs)
    assert res < mp.mpf(10)**-80, res
    cc, E = [sol[j] for j in range(n + 1)], sol[n + 1]
    E = abs(E)
    err = [fx[k] - sum(cc[j]*xs[k]**j for j in range(n + 1)) for k in range(m)]
    mxin = max(abs(err[k]) for k in range(m) if k not in act)
    duals = r.ineqlin.marginals
    return E, act, mxin/E, min(abs(d) for d in duals if abs(d) > 1e-12)


def solve_complex(f, zs, n):
    """Real-symmetric data (f(conj z) = conj f(z), conjugation-closed samples):
    the unique best approximation has REAL coefficients, so the KKT system is
    set up in the n+1 real coefficients over the active set reduced modulo
    conjugation (a pair counted once), which makes it square and regular."""
    fz = [f(w) for w in zs]
    Z = np.array([complex(w) for w in zs]); Fz = np.array([complex(v) for v in fz])
    V = np.vander(Z, n + 1, increasing=True)
    M = len(Z); N2 = 2*(n + 1)
    err = lambda x: Fz - V @ (x[:n + 1] + 1j*x[n + 1:N2])
    x0 = np.zeros(N2 + 1)
    ls = np.linalg.lstsq(V, Fz, rcond=None)[0]; x0[:n + 1] = ls.real; x0[n + 1:N2] = ls.imag
    x0[N2] = np.max(np.abs(err(x0))) * 1.01
    for _ in range(3):
        r = minimize(lambda x: x[N2], x0, constraints=[{"type": "ineq", "fun": lambda x: x[N2]**2 - np.abs(err(x))**2}],
                     method="SLSQP", options={"ftol": 1e-16, "maxiter": 5000})
        x0 = r.x
    E0 = r.x[N2]; e = np.abs(err(r.x))
    act = [k for k in range(M) if e[k] > E0*(1 - 1e-6)]
    red = sorted(set(min(k, (M - k) % M) for k in act))
    m = len(red)
    pw = [[zs[k]**j for j in range(n + 1)] for k in range(M)]

    def ek(a, k):
        return fz[k] - sum(a[j]*pw[k][j] for j in range(n + 1))

    def grads(a, k):
        e_ = ek(a, k)
        return [-2*mp.re(mp.conj(e_)*pw[k][j]) for j in range(n + 1)]

    G = np.array([[float(g) for g in grads([F(v) for v in r.x[:n + 1]], k)] for k in red]).T
    Aug = np.vstack([G, np.ones((1, m))]); rhs = np.r_[np.zeros(n + 1), 1]
    lam0 = np.linalg.lstsq(Aug, rhs, rcond=None)[0]

    def Fsys(*u):
        a, E, lam = u[:n + 1], u[n + 1], u[n + 2:]
        gs = [grads(a, k) for k in red]
        return ([abs(ek(a, k))**2 - E**2 for k in red] +
                [sum(lam[i]*gs[i][j] for i in range(m)) for j in range(n + 1)] + [sum(lam) - 1])

    u0 = [F(v) for v in r.x[:n + 1]] + [F(E0)] + [F(v) for v in lam0]
    u = mp.findroot(Fsys, u0, tol=mp.mpf(10)**-85, maxsteps=100)
    a, E, lam = u[:n + 1], u[n + 1], u[n + 2:]
    redset = set(red)
    mxin = max((abs(ek(a, k)) for k in range(M) if min(k, (M - k) % M) not in redset), default=F(0))
    print("   coeffs", [mp.nstr(v, 20) for v in a], "active", act)
    return abs(E), red, mxin/abs(E), min(lam)


if __name__ == "__main__":
    import sys
    ks = [int(a) for a in sys.argv[1:]] or list(range(1, 11))
    RP, CP = real_probs(), complex_probs()
    for k in ks:
        if k in RP:
            E, act, ratio, dmin = solve_real(*RP[k])
            note = "active %s, max inactive/E %s, min LP dual %.3g" % (act, mp.nstr(ratio, 8), dmin)
        else:
            E, act, ratio, lmin = solve_complex(*CP[k])
            note = "active %s, max inactive/E %s, min lambda %s" % (act, mp.nstr(ratio, 8), mp.nstr(lmin, 5))
        if k == 1:
            note += "; vs 1/8: %s" % mp.nstr(E - F(1)/8, 3)
        print(k, mp.nstr(E, 70), note, flush=True)
