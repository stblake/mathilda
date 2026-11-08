"""Holdout B, family C6 (experiment 231): references for bc6.1 .. bc6.10.

All ten are lattice Green's function problems: G = expected number of visits
(to the origin, or to a site x), return probability 1 - 1/G, hitting
probability G(x)/G(0).  Two independent methods per problem:

  Method A ("I0"), for walks whose steps are +-e_i along coordinate axes
  (product walks, any dimension): Sum_n P_n(x) t^n/n! = Prod_i (p_i/q_i)^(x_i/2)
  I_{x_i}(2 t sqrt(p_i q_i)), so G(x) = Int_0^oo e^-t Prod_i (...) dt (Laplace
  transform of the continuous-time walk; no Fourier integral involved).

  Method B ("1-D"), for any 2-D walk: the Fourier integral with ONE variable
  done in closed form by residues: Int dk/(2 pi) 1/(A - p e^{ik} - q e^{-ik})
  = -1/(p (w_in - w_out)) with w_in the root of p w^2 - A w + q inside the unit
  circle; the other variable by mpmath.quad (tanh-sinh, split at the
  singularities, with k = u^2 and 90 extra digits inside the integrand,
  where A^2 - 4pq cancels near k = 0).  (Check, product walks: the series
  G = sum_n C(2n,n) S_n(EW, NS), n S_n = (2n-1)(x+y) S_{n-1} - (n-1)(x-y)^2
  S_{n-2}, agrees with Method A to 80 digits at e = 0.028.)  For the non-product walks (triangular, king) Method A is
  replaced by Method B with the OTHER variable done in closed form.

  3-D: bc6.5 Glasser-Zucker Gamma closed form (Watson's integral) vs Method A;
  bc6.8 Method A vs a 1-D integral of the complex 2-D square-lattice Green's
  function (2/(pi A)) K(m = 4/(9A^2)) (mpmath.ellipk).

Run: python3.11 bc6_ref.py   (dps 75; ~ a few minutes)
"""
import mpmath as mp
from mpmath import mpf, pi, quad, inf, findroot, besseli, exp, sqrt, cos, sin, nstr

import sys
SKIP = set(sys.argv[1:])   # problem numbers to skip
mp.mp.dps = 75
# step probabilities are made at 400 digits: at the singular point the
# integrand needs 1 - sum(p) = 0 to far more than dps digits, otherwise the
# k = 0 singularity is regularised at scale delta and the integral moves by
# O(sqrt(delta)) (an error of 1e-38 at dps 75).
with mp.workdps(400):
    Q = mpf(1)/4; S = mpf(1)/6; T3 = mpf(1)/3; E8 = mpf(1)/8

def GI(axes, site=None):
    """Method A: Int_0^oo e^-t Prod I_x(2 t sqrt(pq)) (p/q)^(x/2) dt."""
    def f(t):
        r = exp(-t)
        for i, (p, q) in enumerate(axes):
            x = 0 if site is None else site[i]
            r *= sqrt(p/q)**x * besseli(x, 2*t*sqrt(p*q))
        return r
    return quad(f, [0, 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, inf])


def inner(A, p, q):
    """(1/2pi) Int_0^2pi dk / (A - p e^{ik} - q e^{-ik}) by residues."""
    if p == 0:
        return 1/A
    d = sqrt(A*A - 4*p*q)
    w1, w2 = (A + d)/(2*p), (A - d)/(2*p)
    if abs(w1) < abs(w2):
        w1, w2 = w2, w1          # w2 inside
    return -1/(p*(w2 - w1))


def G1D(coef, site1=0, pts=None):
    """Method B: (1/2pi) Int_{-pi}^{pi} dk Re[e^{-i k site1} inner(A(k),p(k),q(k))],
    folded onto [0,pi] and with k = u^2 (removes the k^(-1/2) singularity at 0)."""
    def f(k):
        with mp.extradps(90):    # A^2 - 4pq cancels catastrophically near k = 0
            A, p, q = coef(k)
            return +mp.re(mp.expj(-k*site1)*inner(A, p, q))
    def g(u):
        k = u*u
        return 2*u*(f(k) + f(-k))
    r = sqrt(pi)
    return quad(g, mp.linspace(0, r, 9))/(2*pi)


def EW(k, E, W):
    return E*mp.expj(k) + W*mp.expj(-k)


def show(name, a, b):
    print(name, nstr(a, 68))
    print(' ' * len(name), nstr(b, 68), ' agree to', nstr(-mp.log10(abs(a - b)/abs(a)), 3), 'digits', flush=True)


def root(fn, x0):
    x = mp.re(findroot(fn, mpf(x0), tol=mpf(10)**-72, verify=False))
    print('   residual', nstr(abs(fn(x)), 3), flush=True)
    return x


# --- Z^2 product walks: Method A vs Method B (inner over k2, N/S as p, q)
def sq(E, W, N, Sx):
    """E, W, N, Sx: functions of nothing, evaluated inside the integrand's
    raised precision (so that E + W + N + Sx = 1 to that precision)."""
    a = lambda: GI([(E(), W()), (N(), Sx())])
    b = lambda: G1D(lambda k: (1 - EW(k, E(), W()), N(), Sx()))
    return a, b

# bc6.1  E/W bias, return probability 0.6  (G = 5/2)
fa = lambda e: sq(lambda: Q+e, lambda: Q-e, lambda: Q, lambda: Q)[0]() - mpf(5)/2
fb = lambda e: sq(lambda: Q+e, lambda: Q-e, lambda: Q, lambda: Q)[1]() - mpf(5)/2
if '1' not in SKIP: show('bc6.1', root(fa, '0.02797'), root(fb, '0.02797'))
# bc6.2  E and N both 1/4+e, W and S 1/4-e; return probability 1/2 (G = 2)
fa = lambda e: sq(lambda: Q+e, lambda: Q-e, lambda: Q+e, lambda: Q-e)[0]() - 2
fb = lambda e: sq(lambda: Q+e, lambda: Q-e, lambda: Q+e, lambda: Q-e)[1]() - 2
if '2' not in SKIP: show('bc6.2', root(fa, '0.04394'), root(fb, '0.04394'))
# bc6.3  N,S 1/6; E 1/3+e, W 1/3-e; return probability 1/2
fa = lambda e: sq(lambda: T3+e, lambda: T3-e, lambda: S, lambda: S)[0]() - 2
fb = lambda e: sq(lambda: T3+e, lambda: T3-e, lambda: S, lambda: S)[1]() - 2
show('bc6.3', root(fa, '0.08003'), root(fb, '0.08003'))


# bc6.4  triangular lattice: steps +-e1, +-e2, +-(e1+e2), 1/6 each, +e1 gets +e, -e1 gets -e
def tri_a(e):   # inner over k2
    return G1D(lambda k: (1 - EW(k, S+e, S-e), (1 + mp.expj(k))/6, (1 + mp.expj(-k))/6))
def tri_b(e):   # inner over k1
    return G1D(lambda k: (1 - cos(k)/3, S + e + mp.expj(k)/6, S - e + mp.expj(-k)/6))
show('bc6.4', root(lambda e: tri_a(e) - 2, '0.04637'), root(lambda e: tri_b(e) - 2, '0.04637'))

# bc6.5  Polya's constant (3-D simple cubic, unbiased): 1 - 1/G
from mpmath import gamma
Gz = sqrt(6)/(32*pi**3)*gamma(mpf(1)/24)*gamma(mpf(5)/24)*gamma(mpf(7)/24)*gamma(mpf(11)/24)
GA = GI([(S, S)]*3)
show('bc6.5', 1 - 1/Gz, 1 - 1/GA)


# bc6.6  king walk (8 neighbours, 1/8 each), E 1/8+e, W 1/8-e
def king_a(e):  # inner over k2: p = q = (1 + 2 cos k1)/8
    return G1D(lambda k: (1 - EW(k, E8+e, E8-e), (1 + 2*cos(k))/8, (1 + 2*cos(k))/8))
def king_b(e):  # inner over k1
    return G1D(lambda k: (1 - cos(k)/4, E8 + e + cos(k)/4, E8 - e + cos(k)/4))
show('bc6.6', root(lambda e: king_a(e) - 2, '0.03288'), root(lambda e: king_b(e) - 2, '0.03288'))

# bc6.7  original model at e = 1/10: return probability
e = mpf(1)/10
a, b = sq(lambda: Q+e, lambda: Q-e, lambda: Q, lambda: Q)
show('bc6.7', 1 - 1/a(), 1 - 1/b())


# bc6.8  3-D: E 1/6+e, W 1/6-e, N,S,U,D 1/6: return probability 1/4 (G = 4/3)
def g3b(e):
    def f(k):
        with mp.extradps(90):
            A = 1 - EW(k, S+e, S-e)
            return +mp.re(2/(pi*A)*mp.ellipk(4/(9*A*A)))
    g = lambda u: 2*u*(f(u*u) + f(-u*u))
    return quad(g, mp.linspace(0, sqrt(pi), 9))/(2*pi)
show('bc6.8', root(lambda e: GI([(S+e, S-e), (S, S), (S, S)]) - mpf(4)/3, '0.07344'),
     root(lambda e: g3b(e) - mpf(4)/3, '0.07344'))

# bc6.9  original model at e = 1/20: probability of ever visiting (1,0) = G((1,0))/G(0)
e = mpf(1)/20
pa = GI([(Q+e, Q-e), (Q, Q)], site=(1, 0))/GI([(Q+e, Q-e), (Q, Q)])
pb = (G1D(lambda k: (1 - EW(k, Q+e, Q-e), Q, Q), site1=1)
      / G1D(lambda k: (1 - EW(k, Q+e, Q-e), Q, Q)))
show('bc6.9', pa, pb)

# bc6.10 diagonal-step walk NE 1/4+e, SW 1/4-e, NW, SE 1/4: return probability 0.4 (G = 5/3).
# The rotation (x,y) -> ((x+y)/2, (x-y)/2) maps it onto the Z^2 walk with E/W bias e.
fa = lambda e: sq(lambda: Q+e, lambda: Q-e, lambda: Q, lambda: Q)[0]() - mpf(5)/3
# Method B on the diagonal walk itself: inner over k2 of 1 - cos k1 cos k2 - 2 i e sin(k1+k2):
# e^{ik2} coefficient p = e^{ik1}(1/4+e) + e^{-ik1}/4, q = e^{-ik1}(1/4-e) + e^{ik1}/4, A = 1.
fb = lambda e: G1D(lambda k: (mpf(1), (Q+e)*mp.expj(k) + Q*mp.expj(-k),
                              (Q-e)*mp.expj(-k) + Q*mp.expj(k)),
                   pts=[-pi, -pi/2, 0, pi/2, pi]) - mpf(5)/3
show('bc6.10', root(fa, '0.10608'), root(fb, '0.10608'))
