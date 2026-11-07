# D3 reference, method B: Thomson problem N = 100.
# Parametrisation as in thomson_100.m: point 0 at the north pole, point 1 at
# phi_1 = pi (3 - sqrt 5) (both held), unknowns theta_1..theta_99,
# phi_2..phi_99 (197).  Start: Fibonacci lattice.
#  1. scipy BFGS with the analytic gradient (double) from the Fibonacci start
#     -> E = 4448.350634331..., and 300 random starts confirm that nothing
#     lower is found (global minimum; Wales' Cambridge Cluster Database value
#     4448.350634331).
#  2. Newton in mpmath at dps 130: analytic gradient, Hessian by central
#     differences of the gradient (h = 1e-55), LU solve; the Hessian at the
#     solution is checked positive definite (a strict local minimum).
import numpy as np, scipy.optimize as so, mpmath as mp, sys, time
N = 100
PH1 = np.pi * (3 - np.sqrt(5))
def unpack(v, lib=np):
    th = [0] + list(v[:N-1]); ph = [0, None] + list(v[N-1:])
    return th, ph
def grad_np(v):
    th = np.concatenate([[0.0], v[:N-1]]); ph = np.concatenate([[0.0, PH1], v[N-1:]])
    X = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
    d = X[:, None, :] - X[None, :, :]; r2 = (d**2).sum(-1); np.fill_diagonal(r2, 1.0)
    r = np.sqrt(r2); w = 1/r**3; np.fill_diagonal(w, 0.0)
    gX = -(w[:, :, None]*d).sum(1)
    dth = np.stack([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), -np.sin(th)], 1)
    dph = np.stack([-np.sin(th)*np.sin(ph), np.sin(th)*np.cos(ph), 0*th], 1)
    iu = np.triu_indices(N, 1)
    return (1/r)[iu].sum(), np.concatenate([(gX*dth).sum(1)[1:], (gX*dph).sum(1)[2:]])
k = np.arange(N)
v0 = np.concatenate([np.arccos(1 - 2*k[1:]/N), (k[2:]*PH1)])
r = so.minimize(grad_np, v0, jac=True, method='BFGS', options={'gtol': 1e-11})
print("scipy BFGS from Fibonacci start:", repr(r.fun), flush=True)
rng = np.random.default_rng(7); best = 1e9
for s in range(300):
    v = np.concatenate([np.arccos(rng.uniform(-1, 1, N-1)), rng.uniform(0, 2*np.pi, N-2)])
    rr = so.minimize(grad_np, v, jac=True, method='BFGS', options={'gtol': 1e-9})
    best = min(best, rr.fun)
print("best of 300 random starts:", repr(best), flush=True)

mp.mp.dps = 130
ph1 = mp.pi * (3 - mp.sqrt(5))
def grad_mp(v):
    th = [mp.mpf(0)] + list(v[:N-1]); ph = [mp.mpf(0), ph1] + list(v[N-1:])
    st = [mp.sin(t) for t in th]; ct = [mp.cos(t) for t in th]
    sp = [mp.sin(p) for p in ph]; cp = [mp.cos(p) for p in ph]
    X = [(st[i]*cp[i], st[i]*sp[i], ct[i]) for i in range(N)]
    g = [[mp.mpf(0)]*3 for _ in range(N)]; E = mp.mpf(0)
    for i in range(N):
        xi = X[i]
        for j in range(i+1, N):
            xj = X[j]
            dx = (xi[0]-xj[0], xi[1]-xj[1], xi[2]-xj[2])
            ir = 1/mp.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2); ir3 = ir**3
            E += ir
            for a in range(3):
                g[i][a] -= dx[a]*ir3; g[j][a] += dx[a]*ir3
    gt = [g[i][0]*ct[i]*cp[i] + g[i][1]*ct[i]*sp[i] - g[i][2]*st[i] for i in range(N)]
    gp = [-g[i][0]*st[i]*sp[i] + g[i][1]*st[i]*cp[i] for i in range(N)]
    return E, gt[1:] + gp[2:]
v = [mp.mpf(float(t)) for t in r.x]
n = len(v); h = mp.mpf(10)**-55
for it in range(6):
    t0 = time.time()
    E, g = grad_mp(v)
    gn = mp.sqrt(sum(x**2 for x in g))
    print("Newton", it, "|grad| =", mp.nstr(gn, 5), " E =", mp.nstr(E, 70), flush=True)
    if gn < mp.mpf(10)**-115: break
    H = mp.matrix(n, n)
    for kk in range(n):
        vp = list(v); vp[kk] += h; vm = list(v); vm[kk] -= h
        gpp = grad_mp(vp)[1]; gmm = grad_mp(vm)[1]
        for i in range(n): H[i, kk] = (gpp[i] - gmm[i])/(2*h)
    dv = mp.lu_solve(H, mp.matrix(g))
    v = [v[i] - dv[i] for i in range(n)]
    print("   step time %.0fs" % (time.time()-t0), flush=True)
Hd = np.array([[float(H[i, j]) for j in range(n)] for i in range(n)])
print("Hessian eigenvalues min/max:", np.linalg.eigvalsh((Hd+Hd.T)/2)[[0, -1]])
print("E_100 =", mp.nstr(E, 80))
