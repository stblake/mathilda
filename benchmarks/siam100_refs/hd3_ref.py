"""Holdout A, family hd3 (experiment 208): minimal point energies.

Problem k: (manifold, N, kernel).  Sphere: the parametrisation of thomson_100.m
(point 0 at the north pole, phi_1 = pi(3 - sqrt 5) held, unknowns
theta_1..theta_{N-1}, phi_2..phi_{N-1}), Fibonacci start.  Circle: angles
t_1..t_{N-1} (t_0 = 0 held), start t_k = 2 pi (k/N)^(3/2).
Kernels: 'c' 1/r, ('r', s) r^-s, 'log' -log r, 'dist' -r.
Steps: (1) scipy BFGS with the analytic gradient from the benchmark start;
(2) random-start BFGS to confirm nothing lower exists (global minimum);
(3) mpmath Newton at dps 120 (analytic gradient, central-difference Hessian,
pseudo-inverse step that drops the gauge null directions).
Writes hd3_start_<k>.wl (the double minimiser) for hd3_ref.wl.
Usage: hd3_ref.py k [nrand]"""
import sys, time
import numpy as np, scipy.optimize as so, mpmath as mp

PROBLEMS = {1: ('S', 5, 'c'), 2: ('S', 12, 'c'), 3: ('S', 8, 'c'), 4: ('S', 32, 'c'),
            5: ('S', 20, ('r', 2)), 6: ('S', 30, 'log'), 7: ('S', 16, 'dist'),
            8: ('C', 30, 'c'), 9: ('S', 50, 'c'), 10: ('S', 72, 'c')}


def kern(kind, r, lib):
    """g(r), g'(r)/r."""
    if kind == 'c':
        return 1 / r, -1 / r ** 3
    if kind == 'log':
        return -lib.log(r), -1 / r ** 2
    if kind == 'dist':
        return -r, -1 / r
    s = kind[1]
    return r ** (-s), -s * r ** (-s - 2)


def points_np(v, man, N):
    if man == 'C':
        t = np.concatenate([[0.0], v])
        X = np.stack([np.cos(t), np.sin(t)], 1)
        return X, [np.stack([-np.sin(t), np.cos(t)], 1)], [1]
    PH1 = np.pi * (3 - np.sqrt(5))
    th = np.concatenate([[0.0], v[:N - 1]]); ph = np.concatenate([[0.0, PH1], v[N - 1:]])
    X = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
    dth = np.stack([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)], 1)
    dph = np.stack([-np.sin(th) * np.sin(ph), np.sin(th) * np.cos(ph), 0 * th], 1)
    return X, [dth, dph], [1, 2]


def fg_np(v, man, N, kind):
    X, dirs, offs = points_np(v, man, N)
    d = X[:, None, :] - X[None, :, :]
    r2 = (d ** 2).sum(-1); np.fill_diagonal(r2, 1.0); r = np.sqrt(r2)
    g, gr = kern(kind, r, np); np.fill_diagonal(g, 0.0); np.fill_diagonal(gr, 0.0)
    E = g.sum() / 2
    gX = (gr[:, :, None] * d).sum(1)
    grad = np.concatenate([(gX * D).sum(1)[o:] for D, o in zip(dirs, offs)])
    return E, grad


def start(man, N):
    k = np.arange(N)
    if man == 'C':
        return 2 * np.pi * (k[1:] / N) ** 1.5
    PH1 = np.pi * (3 - np.sqrt(5))
    return np.concatenate([np.arccos(1 - 2 * k[1:] / N), k[2:] * PH1])


def rand_start(man, N, rng):
    if man == 'C':
        return rng.uniform(0, 2 * np.pi, N - 1)
    return np.concatenate([np.arccos(rng.uniform(-1, 1, N - 1)), rng.uniform(0, 2 * np.pi, N - 2)])


def fg_mp(v, man, N, kind):
    if man == 'C':
        t = [mp.mpf(0)] + list(v)
        X = [(mp.cos(a), mp.sin(a)) for a in t]
        dirs = [[(-mp.sin(a), mp.cos(a)) for a in t]]; offs = [1]
    else:
        ph1 = mp.pi * (3 - mp.sqrt(5))
        th = [mp.mpf(0)] + list(v[:N - 1]); ph = [mp.mpf(0), ph1] + list(v[N - 1:])
        st = [mp.sin(a) for a in th]; ct = [mp.cos(a) for a in th]
        sp = [mp.sin(a) for a in ph]; cp = [mp.cos(a) for a in ph]
        X = [(st[i] * cp[i], st[i] * sp[i], ct[i]) for i in range(N)]
        dirs = [[(ct[i] * cp[i], ct[i] * sp[i], -st[i]) for i in range(N)],
                [(-st[i] * sp[i], st[i] * cp[i], mp.mpf(0)) for i in range(N)]]
        offs = [1, 2]
    dim = len(X[0])
    G = [[mp.mpf(0)] * dim for _ in range(N)]; E = mp.mpf(0)
    for i in range(N):
        for j in range(i + 1, N):
            dx = [X[i][a] - X[j][a] for a in range(dim)]
            r = mp.sqrt(sum(x * x for x in dx))
            g, gr = kern(kind, r, mp)
            E += g
            for a in range(dim):
                G[i][a] += gr * dx[a]; G[j][a] -= gr * dx[a]
    grad = []
    for D, o in zip(dirs, offs):
        grad += [sum(G[i][a] * D[i][a] for a in range(dim)) for i in range(o, N)]
    return E, grad


if __name__ == "__main__":
    k = int(sys.argv[1]); nrand = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    man, N, kind = PROBLEMS[k]
    f = lambda v: fg_np(v, man, N, kind)
    r = so.minimize(f, start(man, N), jac=True, method='BFGS', options={'gtol': 1e-11})
    print("BFGS from benchmark start:", repr(r.fun), flush=True)
    r2 = so.minimize(lambda v: f(v)[0], start(man, N), method='BFGS')
    print("scipy minimize default (FD gradient) from start:", repr(r2.fun), flush=True)
    rng = np.random.default_rng(11); vals = []
    for _ in range(nrand):
        rr = so.minimize(f, rand_start(man, N, rng), jac=True, method='BFGS', options={'gtol': 1e-10})
        vals.append(rr.fun)
    vals = np.array(vals)
    print("random starts: best %.15g, frac within 1e-8 of start-min: %.2f" %
          (vals.min(), np.mean(np.abs(vals - r.fun) < 1e-8 * abs(r.fun))), flush=True)
    with open('hd3_start_%d.wl' % k, 'w') as fh:
        fh.write('hstart = {' + ', '.join(repr(float(t)).replace('e', '*^') for t in r.x) + '};\n')
    mp.mp.dps = 120
    v = [mp.mpf(float(t)) for t in r.x]; n = len(v); h = mp.mpf(10) ** -50
    for it in range(8):
        E, g = fg_mp(v, man, N, kind)
        gn = mp.sqrt(sum(x ** 2 for x in g))
        print("Newton", it, "|grad| =", mp.nstr(gn, 5), "E =", mp.nstr(E, 70), flush=True)
        if gn < mp.mpf(10) ** -105:
            break
        H = mp.matrix(n, n)
        for kk in range(n):
            vp = list(v); vp[kk] += h; vm = list(v); vm[kk] -= h
            gp = fg_mp(vp, man, N, kind)[1]; gm = fg_mp(vm, man, N, kind)[1]
            for i in range(n):
                H[i, kk] = (gp[i] - gm[i]) / (2 * h)
        # pseudo-inverse step: a point other than 0 at a pole makes its phi a
        # flat (gauge) direction, so drop Hessian eigenvalues below 1e-5 (the
        # physical ones are >= 4e-4 in every problem here)
        Hs = (H + H.T) / 2
        ew, Q = mp.eigsy(Hs)
        gv = mp.matrix(g)
        dv = mp.matrix(n, 1)
        for i in range(n):
            if abs(ew[i]) > mp.mpf(10) ** -5:
                qi = Q[:, i]
                dv += qi * ((qi.T * gv)[0] / ew[i])
        v = [v[i] - dv[i] for i in range(n)]
    Hd = np.array([[float(H[i, j]) for j in range(n)] for i in range(n)])
    ev = np.linalg.eigvalsh((Hd + Hd.T) / 2)
    print("Hessian eig min/max:", ev[0], ev[-1])
    print("RESULT", k, mp.nstr(E, 80))
