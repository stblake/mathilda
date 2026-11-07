"""Holdout B, family D3 reference, method A (scipy + mpmath).

Sphere problems (bd3.1, 2, 4, 5, 6, 7, 9, 10): N points on the unit sphere,
pair kernel K(r) = 1/r (Coulomb), -log r (log), 1/r^2 (Riesz s=2), -r (sum of
distances, maximised = minimise -sum).  Parametrisation as in the .m file:
point 0 at the north pole, point 1 at phi_1 = pi (3 - sqrt 5), unknowns
theta_1..theta_{N-1}, phi_2..phi_{N-1}; start = Fibonacci lattice.
  1. scipy BFGS with analytic gradient from the Fibonacci start (double), and
     NRAND random starts (BFGS, gtol 1e-10): nothing lower may be found.
  2. Newton in mpmath (dps 80): analytic gradient, Hessian by central
     differences (h = 1e-35), LU; positive-definite Hessian checked.
bd3.3 (30 points on the unit circle, Coulomb; angles phi_1..phi_29, phi_0 = 0,
start phi_k = 2 pi (k/N)^2) and bd3.8 (Lennard-Jones 13, 39 Cartesian
coordinates, start = centre + 12 Fibonacci points at radius 1.1) use the same
two steps (LJ13: Newton on the 1-parameter icosahedral reduction
E(r) = 12 V(r) + 30 V(a) + 30 V(phi a) + 6 V(2 r), a = edge = r / sin(2 pi/5),
after checking the BFGS minimiser IS the icosahedron).
Writes the double minimiser to bd3_<p>.start for method B (bd3_ref.gp).
Usage: bd3_ref.py p [NRAND]"""
import sys, time
import numpy as np, scipy.optimize as so, mpmath as mp

p = int(sys.argv[1]); NRAND = int(sys.argv[2]) if len(sys.argv) > 2 else 30
SPEC = {1: (5, 'coul'), 2: (12, 'coul'), 4: (24, 'coul'), 5: (32, 'coul'), 6: (20, 'log'),
        7: (20, 'r2'), 9: (60, 'coul'), 10: (30, 'dist')}
PH = np.pi * (3 - np.sqrt(5))
# K(r), K'(r)
KN = {'coul': (lambda r: 1/r, lambda r: -1/r**2), 'log': (lambda r: -np.log(r), lambda r: -1/r),
      'r2': (lambda r: 1/r**2, lambda r: -2/r**3), 'dist': (lambda r: -r, lambda r: -1 + 0*r)}
KM = {'coul': (lambda r: 1/r, lambda r: -1/r**2), 'log': (lambda r: -mp.log(r), lambda r: -1/r),
      'r2': (lambda r: 1/r**2, lambda r: -2/r**3), 'dist': (lambda r: -r, lambda r: mp.mpf(-1))}


def newton(grad, v, n, tol, label, h=mp.mpf(10)**-35, free=None):
    # free: the coordinates Newton moves; a phi of a point sitting at the
    # south pole is a coordinate singularity (zero Hessian mode) and is frozen.
    free = list(range(n)) if free is None else free
    nf = len(free); n_all = n; n = nf
    for it in range(12):
        t0 = time.time()
        E, g = grad(v)
        g = [g[i] for i in free]
        gn = mp.sqrt(sum(x**2 for x in g))
        print("Newton", it, "|grad| =", mp.nstr(gn, 5), " E =", mp.nstr(E, 70), flush=True)
        if gn < tol:
            break
        H = mp.matrix(n, n)
        for kk in range(n):
            k = free[kk]
            vp = list(v); vp[k] += h; vm = list(v); vm[k] -= h
            gpp = grad(vp)[1]; gmm = grad(vm)[1]
            for ii in range(n):
                H[ii, kk] = (gpp[free[ii]] - gmm[free[ii]]) / (2*h)
        dv = mp.lu_solve(H, mp.matrix(g))
        v = list(v)
        for ii in range(n):
            v[free[ii]] -= dv[ii]
        print("   step %.0fs" % (time.time() - t0), flush=True)
    Hd = np.array([[float(H[i, j]) for j in range(n)] for i in range(n)])
    ev = np.linalg.eigvalsh((Hd + Hd.T)/2)
    print("Hessian eigenvalues min/max:", ev[0], ev[-1])
    print("%s %s" % (label, mp.nstr(E, 70)))
    return v


mp.mp.dps = 80
if p in SPEC:
    N, kind = SPEC[p]; K, dK = KN[kind]; Km, dKm = KM[kind]

    def fg(v):
        th = np.concatenate([[0.0], v[:N-1]]); ph = np.concatenate([[0.0, PH], v[N-1:]])
        X = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
        d = X[:, None, :] - X[None, :, :]; r2 = (d**2).sum(-1); np.fill_diagonal(r2, 1.0)
        r = np.sqrt(r2); w = dK(r)/r; np.fill_diagonal(w, 0.0)
        gX = (w[:, :, None]*d).sum(1)
        dth = np.stack([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), -np.sin(th)], 1)
        dph = np.stack([-np.sin(th)*np.sin(ph), np.sin(th)*np.cos(ph), 0*th], 1)
        iu = np.triu_indices(N, 1)
        return K(r)[iu].sum(), np.concatenate([(gX*dth).sum(1)[1:], (gX*dph).sum(1)[2:]])
    k = np.arange(N)
    v0 = np.concatenate([np.arccos(1 - 2*k[1:]/N), k[2:]*PH])
    r = so.minimize(fg, v0, jac=True, method='BFGS', options={'gtol': 1e-11})
    print("BFGS from Fibonacci start:", repr(r.fun), flush=True)
    rng = np.random.default_rng(11); vals = []
    for s in range(NRAND):
        v = np.concatenate([np.arccos(rng.uniform(-1, 1, N-1)), rng.uniform(0, 2*np.pi, N-2)])
        vals.append(so.minimize(fg, v, jac=True, method='BFGS', options={'gtol': 1e-10}).fun)
    vals = np.sort(vals)
    print("random starts: best %.13f, fraction within 1e-8 of Fibonacci-start value %.2f"
          % (vals[0], np.mean(np.abs(vals - r.fun) < 1e-8)), flush=True)
    np.savetxt("bd3_%d.start" % p, r.x, fmt="%.17g")
    ph1 = mp.pi * (3 - mp.sqrt(5))

    def gm(v):
        th = [mp.mpf(0)] + list(v[:N-1]); ph = [mp.mpf(0), ph1] + list(v[N-1:])
        st = [mp.sin(t) for t in th]; ct = [mp.cos(t) for t in th]
        sp = [mp.sin(q) for q in ph]; cp = [mp.cos(q) for q in ph]
        X = [(st[i]*cp[i], st[i]*sp[i], ct[i]) for i in range(N)]
        g = [[mp.mpf(0)]*3 for _ in range(N)]; E = mp.mpf(0)
        for i in range(N):
            for j in range(i+1, N):
                dx = [X[i][a] - X[j][a] for a in range(3)]
                rr = mp.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)
                E += Km(rr); w = dKm(rr)/rr
                for a in range(3):
                    g[i][a] += w*dx[a]; g[j][a] -= w*dx[a]
        gt = [g[i][0]*ct[i]*cp[i] + g[i][1]*ct[i]*sp[i] - g[i][2]*st[i] for i in range(N)]
        gq = [-g[i][0]*st[i]*sp[i] + g[i][1]*st[i]*cp[i] for i in range(N)]
        return E, gt[1:] + gq[2:]
    frozen = [N - 1 + (i - 2) for i in range(2, N) if abs(np.sin(r.x[i - 1])) < 1e-4]
    print("frozen phi coordinates (points at the south pole):", frozen)
    newton(gm, [mp.mpf(float(x)) for x in r.x], 2*N - 3, mp.mpf(10)**-70, "bd3.%d" % p,
           free=[i for i in range(2*N - 3) if i not in frozen])
elif p == 3:
    N = 30

    def fg(v):
        ph = np.concatenate([[0.0], v]); d = ph[:, None] - ph[None, :]
        r = 2*np.abs(np.sin(d/2)); np.fill_diagonal(r, 1.0)
        dr = np.sign(np.sin(d/2))*np.cos(d/2); w = -dr/r**2; np.fill_diagonal(w, 0.0)
        iu = np.triu_indices(N, 1)
        return (1/r)[iu].sum(), w.sum(1)[1:]
    k = np.arange(1, N)
    r = so.minimize(fg, 2*np.pi*(k/N)**2, jac=True, method='BFGS', options={'gtol': 1e-11})
    print("BFGS:", repr(r.fun))

    def gm(v):
        ph = [mp.mpf(0)] + list(v); E = mp.mpf(0); g = [mp.mpf(0)]*N
        for i in range(N):
            for j in range(i+1, N):
                d = ph[i] - ph[j]; rr = 2*mp.sin(d/2)
                if rr < 0: rr = -rr; s = -1
                else: s = 1
                E += 1/rr; w = -s*mp.cos(d/2)/rr**2
                g[i] += w; g[j] -= w
        return E, g[1:]
    newton(gm, [mp.mpf(float(x)) for x in r.x], N-1, mp.mpf(10)**-70, "bd3.3")
    print("closed form (N/2) sum 1/(2 sin(pi k/N)):",
          mp.nstr(mp.mpf(N)/2*mp.fsum(1/(2*mp.sin(mp.pi*kk/N)) for kk in range(1, N)), 70))
elif p == 8:
    def fg(x):
        X = x.reshape(13, 3); d = X[:, None, :] - X[None, :, :]; r2 = (d**2).sum(-1)
        np.fill_diagonal(r2, 1.0); i6 = 1/r2**3; V = 4*(i6**2 - i6); np.fill_diagonal(V, 0)
        dV = 4*(-12*i6**2 + 6*i6)/r2; np.fill_diagonal(dV, 0)
        return V.sum()/2, (dV[:, :, None]*d).sum(1).ravel()
    k = np.arange(12); th = np.arccos(1 - (2*k + 1)/12); ph = k*PH
    X0 = np.vstack([[0, 0, 0], 1.1*np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)])
    r = so.minimize(fg, X0.ravel(), jac=True, method='BFGS', options={'gtol': 1e-11})
    X = r.x.reshape(13, 3); c = X[np.argmin(np.linalg.norm(X - X.mean(0), axis=1))]
    rad = np.sort(np.linalg.norm(X - c, axis=1))[1:]
    print("BFGS:", repr(r.fun), " radii spread from centre atom:", rad.min(), rad.max())
    np.savetxt("bd3_8.start", r.x, fmt="%.17g")
    phi = (1 + mp.sqrt(5))/2
    V = lambda s: 4*(s**-12 - s**-6)
    E = lambda rr: 12*V(rr) + 30*V(rr/mp.sin(2*mp.pi/5)) + 30*V(phi*rr/mp.sin(2*mp.pi/5)) + 6*V(2*rr)
    rr = mp.findroot(lambda t: mp.diff(E, t), mp.mpf(float(rad.mean())))
    print("icosahedral radius r* =", mp.nstr(rr, 40), " E''(r*) =", mp.nstr(mp.diff(E, rr, 2), 8))
    print("bd3.8", mp.nstr(E(rr), 70))
