"""hc4 (holdout A, family C4) reference, method A: exhaustive grid + Newton.
For each objective: a dense numpy grid over the box (2-D: 4001^2 points, plus
the four box edges at 40001 points; 3-D separable k = 3 done per coordinate),
the 40 best grid cells polished by scipy L-BFGS-B, the best few then by
mpmath.findroot on the gradient (mpmath.diff) at dps 90 -- or, for a
minimiser on the box boundary (k = 8, x = 512), on d/dy with x fixed (and
df/dx < 0 checked).  Closed forms for k = 2 (5/(4 pi)) and k = 3 (3 g(t*),
4 t*^3 - 32 t* + 5 = 0).  Method B: hc4_ref.wl (NMinimize DifferentialEvolution
+ FindMinimum at WorkingPrecision 120)."""
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

mp.mp.dps = 90


def fns(L):
    """L = numpy or mpmath module-like namespace"""
    s, c, e, sq, ab, pi = L.sin, L.cos, L.exp, L.sqrt, L.fabs if L is mp else np.abs, L.pi
    return {
        1: (lambda x, y: (4 - 2.1*x**2 + x**4/3)*x**2 + x*y + (-4 + 4*y**2)*y**2
            if L is np else (4 - mp.mpf(21)/10*x**2 + x**4/3)*x**2 + x*y + (-4 + 4*y**2)*y**2,
            [(-3, 3), (-2, 2)]),
        2: (lambda x, y: (y - 5.1/(4*pi**2)*x**2 + 5/pi*x - 6)**2 + 10*(1 - 1/(8*pi))*c(x) + 10
            if L is np else (y - mp.mpf(51)/10/(4*pi**2)*x**2 + 5/pi*x - 6)**2 + 10*(1 - 1/(8*pi))*c(x) + 10,
            [(-5, 10), (0, 15)]),
        4: (lambda x, y: sum(i*c((i + 1)*x + i) for i in range(1, 6))*sum(i*c((i + 1)*y + i) for i in range(1, 6)),
            [(-10, 10), (-10, 10)]),
        5: (lambda x, y: -ab(s(x)*c(y)*e(ab(1 - sq(x*x + y*y)/pi))), [(-10, 10), (-10, 10)]),
        6: (lambda x, y: -(ab(s(x)*s(y)*e(ab(100 - sq(x*x + y*y)/pi))) + 1)**0.1/10000
            if L is np else -(ab(s(x)*s(y)*e(ab(100 - sq(x*x + y*y)/pi))) + 1)**(mp.mpf(1)/10)/10000,
            [(-10, 10), (-10, 10)]),
        7: (lambda x, y: e(s(30*x)) + s(40*e(y)) + s(50*s(x)) + s(s(60*y)) - s(10*(x + y)) + (x*x + y*y)/4,
            [(-1, 1), (-1, 1)]),
        8: (lambda x, y: -(y + 47)*s(sq(ab(x/2 + y + 47))) - x*s(sq(ab(x - (y + 47)))),
            [(-512, 512), (-512, 512)]),
        9: (lambda x, y: (x*x + y - 11)**2 + (x + y*y - 7)**2 + x + y, [(-5, 5), (-5, 5)]),
        10: (lambda x, y: e(s(70*x)) + s(80*e(y)) + s(90*s(x)) + s(s(100*y)) - s(12*(x + y)) + (x*x + y*y)/4,
             [(-1, 1), (-1, 1)]),
    }


NP, MP = fns(np), fns(mp)


def solve2(k, N=4001):
    f, B = NP[k]
    (a, b), (c, d) = B
    X, Y = np.meshgrid(np.linspace(a, b, N), np.linspace(c, d, N), indexing="ij")
    Z = f(X, Y)
    idx = np.argsort(Z, axis=None)[:40]
    pts = [(X.flat[i], Y.flat[i]) for i in idx]
    # the box edges
    t = np.linspace(0, 1, 40001)
    for ex in ([a + 0*t, c + (d - c)*t], [b + 0*t, c + (d - c)*t], [a + (b - a)*t, c + 0*t], [a + (b - a)*t, d + 0*t]):
        z = f(ex[0], ex[1]); i = int(np.argmin(z)); pts.append((ex[0][i], ex[1][i]))
    best = []
    for p in pts:
        r = minimize(lambda v: f(v[0], v[1]), p, method="L-BFGS-B", bounds=B, options={"ftol": 1e-15, "gtol": 1e-12})
        best.append((r.fun, tuple(r.x)))
    best.sort()
    fm = MP[k][0]
    out = []
    for val, (x0, y0) in best[:6]:
        onb = [abs(x0 - a) < 1e-9, abs(x0 - b) < 1e-9, abs(y0 - c) < 1e-9, abs(y0 - d) < 1e-9]
        if any(onb):
            if onb[0] or onb[1]:
                X0 = mp.mpf(a if onb[0] else b)
                yy = mp.findroot(lambda y: mp.diff(lambda u: fm(X0, u), y), mp.mpf(y0))
                xx = X0
                gx = mp.diff(lambda u: fm(u, yy), xx)
                note = "boundary x=%s, df/dx=%s" % (X0, mp.nstr(gx, 5))
            else:
                Y0 = mp.mpf(c if onb[2] else d)
                xx = mp.findroot(lambda x: mp.diff(lambda u: fm(u, Y0), x), mp.mpf(x0))
                yy = Y0
                note = "boundary y"
        else:
            g = [lambda x, y: mp.diff(fm, (x, y), (1, 0)), lambda x, y: mp.diff(fm, (x, y), (0, 1))]
            xx, yy = mp.findroot(g, (mp.mpf(x0), mp.mpf(y0)))
            note = "interior"
        out.append((fm(xx, yy), xx, yy, note, val))
    out.sort(key=lambda r: r[0])
    return out


if __name__ == "__main__":
    import sys
    ks = [int(a) for a in sys.argv[1:]] or [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for k in ks:
        if k == 3:
            # Styblinski-Tang 3-D, separable: 3 * min g, g(t) = (t^4 - 16 t^2 + 5 t)/2 on [-5, 5]
            g = lambda t: (t**4 - 16*t**2 + 5*t)/2
            ts = np.linspace(-5, 5, 400001); t0 = ts[np.argmin(g(ts))]
            tt = mp.findroot(lambda t: 4*t**3 - 32*t + 5, mp.mpf(t0))
            print(3, mp.nstr(3*g(tt), 70), "t* =", mp.nstr(tt, 25), flush=True)
            continue
        out = solve2(k)
        v, x, y, note, v0 = out[0]
        extra = ""
        if k == 2:
            extra = " closed 5/(4pi) agrees %.1f" % float(-mp.log10(abs(v - 5/(4*mp.pi))/v))
        print(k, mp.nstr(v, 70), "at", mp.nstr(x, 20), mp.nstr(y, 20), note,
              "| next:", mp.nstr(out[1][0], 20) if len(out) > 1 else "", extra, flush=True)
