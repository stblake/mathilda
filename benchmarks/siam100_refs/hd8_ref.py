# Holdout-A D8 family (experiment 213) references, method A (mpmath, dps 80).
#   c^2 = (2 pi)^-N Int_{R^N} 1/P(x) dx,  c = sqrt(c^2).
# Sobolev rows P = Sum_{|k|<=m} x^(2k): spherical reduction, radial integral
# in CLOSED FORM by residues (t = r^2, Q(t) = Sum_j h_j(u^2) t^j along the
# direction u), angles by tanh-sinh over one orthant (the integrand is even in
# every coordinate):
#   N = 2:  Int_0^inf r/Q(r^2) dr = (1/2) Int_0^inf dt/Q = -(1/2) Sum_i log(-t_i)/Q'(t_i)
#   N = 3:  Int_0^inf r^2/Q(r^2) dr = (1/2) Int_0^inf t^(1/2)/Q dt
#                                   = (1/2) pi/sin(3 pi/2) Sum_i (-t_i)^(1/2)/Q'(t_i)
# (t_i the roots of Q, principal branches); the formulas are checked against a
# numerical radial quad at a few directions.  1-D rows: tanh-sinh over R.
# Bessel-potential / anisotropic rows: Gamma closed forms vs numerical quad.
# Method B: hd8_ref.wl.
import mpmath as mp
from itertools import product
mp.mp.dps = 80
pi = mp.pi

def hcomplete(y, m):
    """h_0..h_m of the variables y (complete homogeneous symmetric)."""
    p = [None] + [sum(v ** j for v in y) for j in range(1, m + 1)]
    h = [mp.mpf(1)] + [0] * m
    for j in range(1, m + 1):
        h[j] = sum(p[i] * h[j - i] for i in range(1, j + 1)) / j
    return h

def radial(u, m, N, check=False):
    h = hcomplete([v * v for v in u], m)       # Q(t) = sum h_j t^j
    coef = list(reversed(h))                   # highest first for polyroots
    roots = mp.polyroots(coef, maxsteps=200, extraprec=200)
    dQ = lambda t: sum(j * h[j] * t ** (j - 1) for j in range(1, m + 1))
    if N == 2:
        val = -sum(mp.log(-t) / dQ(t) for t in roots) / 2
    else:
        val = (pi / mp.sin(3 * pi / 2)) * sum(mp.sqrt(-t) / dQ(t) for t in roots) / 2
    val = mp.re(val)
    if check:
        Q = lambda t: sum(h[j] * t ** j for j in range(m + 1))
        num = mp.quad(lambda r: r ** (N - 1) / Q(r * r), [0, 1, mp.inf])
        print("     radial check", mp.nstr(val - num, 3))
    return val

def sobolev(m, N):
    if N == 2:
        for th in (mp.mpf('0.3'), mp.mpf('1.1')):
            radial([mp.cos(th), mp.sin(th)], m, 2, check=True)
        I = 4 * mp.quad(lambda th: radial([mp.cos(th), mp.sin(th)], m, 2), [0, pi / 4, pi / 2])
    else:
        radial([mp.mpf('0.6'), mp.mpf('0.0'), mp.mpf('0.8')], m, 3, check=True)
        def F(th, ph):
            st = mp.sin(th)
            return radial([st * mp.cos(ph), st * mp.sin(ph), mp.cos(th)], m, 3) * st
        I = 8 * mp.quad(F, [0, pi / 2], [0, pi / 2])
    return mp.sqrt(I / (2 * pi) ** N)

def c1d(m):
    P = lambda x: sum(x ** (2 * k) for k in range(m + 1))
    return mp.sqrt(mp.quad(lambda x: 1 / P(x), [-mp.inf, -1, 0, 1, mp.inf]) / (2 * pi))

R = {}
R[1] = c1d(1); R[2] = c1d(2); R[3] = c1d(5); R[4] = c1d(10)
print("1-D done", flush=True)
R[5] = sobolev(3, 2); print("5", mp.nstr(R[5], 70), flush=True)
R[6] = sobolev(4, 2); print("6", mp.nstr(R[6], 70), flush=True)
mp.mp.dps = 64      # the 3-D angular quad is the expensive one
R[7] = sobolev(2, 3); print("7", mp.nstr(R[7], 60), flush=True)
mp.mp.dps = 80
# Bessel potential (1 + |x|^2)^-s: Int = pi^(N/2) Gamma(s - N/2)/Gamma(s); numerically by the radial quad
bes = lambda N, s: mp.sqrt(pi ** (mp.mpf(N) / 2) * mp.gamma(s - mp.mpf(N) / 2) / mp.gamma(s) / (2 * pi) ** N)
besnum = lambda N, s, area: mp.sqrt(area * mp.quad(lambda r: r ** (N - 1) * (1 + r * r) ** (-s), [0, 1, mp.inf]) / (2 * pi) ** N)
R[8] = bes(2, mp.mpf(3) / 2); print("8 num-closed", mp.nstr(R[8] - besnum(2, mp.mpf(3) / 2, 2 * pi), 3))
R[9] = bes(3, mp.mpf(5) / 2); print("9 num-closed", mp.nstr(R[9] - besnum(3, mp.mpf(5) / 2, 4 * pi), 3), " 1/(pi sqrt 6):", mp.nstr(R[9] - 1 / (pi * mp.sqrt(6)), 3))
# anisotropic 1/(1 + x^2 + y^4): x-integral pi/sqrt(1+y^4); Int_R dy/sqrt(1+y^4) = Gamma(1/4)^2/(2 sqrt(pi))
R[10] = mp.sqrt(pi * mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(pi)) / (2 * pi) ** 2)
num10 = mp.sqrt(mp.quad(lambda y: pi / mp.sqrt(1 + y ** 4), [-mp.inf, 0, mp.inf]) / (2 * pi) ** 2)
print("10 num-closed", mp.nstr(R[10] - num10, 3))
print("closed 1:", mp.nstr(R[1] - 1 / mp.sqrt(2), 3), " 2:", mp.nstr(R[2] - mp.sqrt(1 / (2 * mp.sqrt(3))), 3))
for k in range(1, 11):
    print("hd8.%d  %s" % (k, mp.nstr(R[k], 70)), flush=True)
