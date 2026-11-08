"""Holdout B, family C8 (experiment 233): references for bc8.1 .. bc8.10.
Each problem: the time t at which a heat-equation solution (separated
eigen-series, as stated in the benchmark) reaches a level.  Two independent
methods per problem:
  A  the eigen-series of the benchmark statement (mpmath.nsum) + findroot;
  B  an independent representation of the same solution:
       image (method-of-images, erfc) series  -- bc8.1, 2, 3, 7, 8, 10
       Jacobi theta transformation             -- bc8.5
       numerical Laplace inversion (mpmath.invertlaplace, Talbot) of the
       exact Laplace-domain solution            -- bc8.4, 6, 9
Run: python3.11 bc8_ref.py   (dps 80)
"""
from mpmath import (mp, mpf, sqrt, pi, exp, erfc, nsum, inf, findroot, nstr, sin, cosh, sinh,
                    besseljzero, besselj, besseli, invertlaplace, log10, tanh)
mp.dps = 80


def show(name, a, b):
    print(name, nstr(a, 70))
    print(' ' * len(name), nstr(b, 70), ' agree to', nstr(-log10(abs(a - b) / abs(a)), 3), 'digits', flush=True)


def root(f, x0):
    x = findroot(f, mpf(x0), tol=mpf(10) ** -78, verify=False)
    print('   residual', nstr(abs(f(x)), 3), flush=True)
    return x


# centre value of the 1-D problem on [-L, L]... here: on an interval of length 2 (half-length 1),
# zero ends, initial 1:  eigen series and image series.
wE = lambda t: 4 / pi * nsum(lambda k: (-1) ** int(k) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi ** 2 * t / 4), [0, inf])
wI = lambda t: 1 - 2 * nsum(lambda n: (-1) ** int(n) * erfc((2 * n + 1) / (2 * sqrt(t))), [0, inf])
# same on an interval of length 1 (half-length 1/2): w1(t) = w(4t)
w1E = lambda t: wE(4 * t)
w1I = lambda t: wI(4 * t)

# bc8.1  square, one side at 5, centre reaches 0.3
show('bc8.1', root(lambda t: mpf(5) / 4 * (1 - wE(t) ** 2) - mpf(3) / 10, '0.13'),
     root(lambda t: mpf(5) / 4 * (1 - wI(t) ** 2) - mpf(3) / 10, '0.13'))

# bc8.2  rod [0,1], u(0)=0, u(1)=1, midpoint reaches 1/4
u2E = lambda t: mpf(1) / 2 - 2 / pi * nsum(lambda k: (-1) ** int(k) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi ** 2 * t), [0, inf])
u2I = lambda t: nsum(lambda n: erfc((2 * n + 1 - mpf(1) / 2) / (2 * sqrt(t))) - erfc((2 * n + 1 + mpf(1) / 2) / (2 * sqrt(t))), [0, inf])
show('bc8.2', root(lambda t: u2E(t) - mpf(1) / 4, '0.05'), root(lambda t: u2I(t) - mpf(1) / 4, '0.05'))

# bc8.3  cube, all faces at 1, centre reaches 9/10
show('bc8.3', root(lambda t: 1 - wE(t) ** 3 - mpf(9) / 10, '0.41'), root(lambda t: 1 - wI(t) ** 3 - mpf(9) / 10, '0.41'))

# bc8.4  unit disk, rim at 1, centre reaches 1/2: u = 1 - 2 sum exp(-j_k^2 t)/(j_k J1(j_k));
# Laplace domain: U(s) = 1/(s I0(sqrt s))
u4E = lambda t: 1 - 2 * nsum(lambda k: exp(-besseljzero(0, int(k)) ** 2 * t) / (besseljzero(0, int(k)) * besselj(1, besseljzero(0, int(k)))), [1, inf])
u4L = lambda t: invertlaplace(lambda s: 1 / (s * besseli(0, sqrt(s))), t, method='talbot')
show('bc8.4', root(lambda t: u4E(t) - mpf(1) / 2, '0.2'), root(lambda t: u4L(t) - mpf(1) / 2, '0.2'))

# bc8.5  unit ball, surface at 1, centre reaches 1/2: u = 1 + 2 sum (-1)^n exp(-n^2 pi^2 t)
#        theta transform: sum_{n in Z} (-1)^n e^{-n^2 pi^2 t} = (pi t)^(-1/2) sum_{n in Z} e^{-(n+1/2)^2/t}
u5E = lambda t: 1 + 2 * nsum(lambda n: (-1) ** int(n) * exp(-n ** 2 * pi ** 2 * t), [1, inf])
u5T = lambda t: 2 / sqrt(pi * t) * nsum(lambda n: exp(-(n + mpf(1) / 2) ** 2 / t), [0, inf])
show('bc8.5', root(lambda t: u5E(t) - mpf(1) / 2, '0.14'), root(lambda t: u5T(t) - mpf(1) / 2, '0.14'))

# bc8.6  square, one side at 5, MEAN temperature reaches 1: (5/4)(1 - wbar^2),
#        wbar = 8/pi^2 sum exp(-(2k+1)^2 pi^2 t/4)/(2k+1)^2 ; Laplace: 1/s - tanh(sqrt s)/s^(3/2)
wbE = lambda t: 8 / pi ** 2 * nsum(lambda k: exp(-(2 * k + 1) ** 2 * pi ** 2 * t / 4) / (2 * k + 1) ** 2, [0, inf])
wbL = lambda t: invertlaplace(lambda s: 1 / s - tanh(sqrt(s)) / s ** (mpf(3) / 2), t, method='talbot')
show('bc8.6', root(lambda t: mpf(5) / 4 * (1 - wbE(t) ** 2) - 1, '0.24'), root(lambda t: mpf(5) / 4 * (1 - wbL(t) ** 2) - 1, '0.24'))

# bc8.7  rectangle [0,2]x[0,1], all sides at 1, centre reaches 1/2: 1 - w(t) w1(t)
show('bc8.7', root(lambda t: 1 - wE(t) * w1E(t) - mpf(1) / 2, '0.09'),
     root(lambda t: 1 - wI(t) * w1I(t) - mpf(1) / 2, '0.09'))

# bc8.8  rod [0,1], ends at 1, centre reaches 1/10, stated as the IMAGE series
#        2 sum (-1)^n erfc((2n+1)/(4 sqrt t)) ; method B: eigen series 1 - w1(t)
u8I = lambda t: 2 * nsum(lambda n: (-1) ** int(n) * erfc((2 * n + 1) / (4 * sqrt(t))), [0, inf])
show('bc8.8', root(lambda t: u8I(t) - mpf(1) / 10, '0.0325'), root(lambda t: 1 - w1E(t) - mpf(1) / 10, '0.0325'))

# bc8.9  square [0,1]^2, all sides at 1, point (1/4,1/4) reaches 1/2: 1 - v^2,
#        v = 4/pi sum sin((2k+1)pi/4)/(2k+1) exp(-(2k+1)^2 pi^2 t); Laplace: 1/s - cosh(sqrt s/4)/(s cosh(sqrt s/2))
vE = lambda t: 4 / pi * nsum(lambda k: sin((2 * k + 1) * pi / 4) / (2 * k + 1) * exp(-(2 * k + 1) ** 2 * pi ** 2 * t), [0, inf])
vL = lambda t: invertlaplace(lambda s: 1 / s - cosh(sqrt(s) / 4) / (s * cosh(sqrt(s) / 2)), t, method='talbot')
show('bc8.9', root(lambda t: 1 - vE(t) ** 2 - mpf(1) / 2, '0.03'), root(lambda t: 1 - vL(t) ** 2 - mpf(1) / 2, '0.03'))

# bc8.10 square, one side at 5, centre reaches 10^-6
show('bc8.10', root(lambda t: mpf(5) / 4 * (1 - wE(t) ** 2) - mpf(10) ** -6, '0.0186'),
     root(lambda t: mpf(5) / 4 * (1 - wI(t) ** 2) - mpf(10) ** -6, '0.0186'))
