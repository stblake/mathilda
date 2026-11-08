"""Holdout B, family C10 (experiment 235): references for bc10.1 .. bc10.10.
Each problem: a Brownian exit (harmonic-measure) probability = the value of a
harmonic function, stated in the benchmark as a separated (Fourier) series.
Two independent methods per problem:
  A  the benchmark's series (mpmath.nsum, or a direct sum to a proven tail
     bound for the double series);
  B  the COMPLEMENTARY expansion: 1 - (probability of the other exits), each
     computed from the separation in the OTHER coordinate (a different,
     slower-converging series, summed at higher precision to absorb the
     cancellation), or a closed form (bc10.5, bc10.7).
Run: python3.11 bc10_ref.py   (dps 120)
"""
from mpmath import (mp, mpf, pi, sech, cosh, sinh, sin, cos, exp, sqrt, atan, nsum, inf, nstr, log10,
                    besseljzero, besselj, besseli, fsum)
mp.dps = 120


def show(name, a, b):
    print(name, nstr(a, 70))
    print(' ' * len(name), nstr(b, 70), ' agree to', nstr(-log10(abs(a - b) / abs(a)), 3), 'digits', flush=True)


sg = lambda k: mpf(-1) if int(k) % 2 else mpf(1)   # mpf: an int quotient would be a float
NS = lambda f: nsum(f, [0, inf])


def ssum(f, tol=mpf(10) ** -125, kmax=200000):
    """direct sum of a series whose terms decay geometrically (monotone envelope)."""
    s = mpf(0); k = 0; small = 0
    while k < kmax:
        t = f(k); s += t
        small = small + 1 if abs(t) < tol else 0
        if small > 5: break
        k += 1
    return s


def rect_ends(L, x0=0):
    """rectangle [-L/2, L/2] x [-1/2, 1/2], ends at 1: value at (x0, 0)."""
    return 4 / pi * NS(lambda k: sg(k) / (2 * k + 1) * cosh((2 * k + 1) * pi * x0) / cosh((2 * k + 1) * pi * L / 2))


def rect_sides(L, x0=0):
    """sides y = +-1/2 at 1: modes cos((2k+1) pi x / L), cosh((2k+1) pi y/L)/cosh((2k+1) pi/(2L))."""
    return 4 / pi * ssum(lambda k: sg(k) / (2 * k + 1) * cos((2 * k + 1) * pi * x0 / L) * sech((2 * k + 1) * pi / (2 * L)))


# bc10.1  4 x 1 rectangle, centre, P(ends)
show('bc10.1', rect_ends(4), 1 - rect_sides(4))
# bc10.2  10 x 1 rectangle, start (3, 0), P(ends)
show('bc10.2', rect_ends(10, 3), 1 - rect_sides(10, 3))


# bc10.3  [0,2] x [0,1], start (1/2,1/2), P(far end x = 2)
def sq_side(W, H, x0, y0):
    """[0,W] x [0,H], side y = H at 1: sum_{n odd} 4/(n pi) sin(n pi x0/W) sinh(n pi y0/W)/sinh(n pi H/W)."""
    return 4 / pi * ssum(lambda k: sin((2 * k + 1) * pi * x0 / W) / (2 * k + 1)
                         * sinh((2 * k + 1) * pi * y0 / W) / sinh((2 * k + 1) * pi * H / W))
A3 = 4 / pi * NS(lambda k: sg(k) / (2 * k + 1) * sinh((2 * k + 1) * pi / 2) / sinh(2 * (2 * k + 1) * pi))
B3 = 1 - sq_side(1, 2, mpf(1) / 2, mpf(3) / 2) - 2 * sq_side(2, 1, mpf(1) / 2, mpf(1) / 2)
show('bc10.3', A3, B3)


# bc10.4  box [-1/2,1/2] x [-1,1] x [-3/2,3/2], centre, P(the two 1 x 2 ends z = +-3/2)
def box_ends(a, b, L):
    """cross-section [-a/2,a/2]x[-b/2,b/2], ends z = +-L/2: double series."""
    c = 16 / pi ** 2
    def row(m):
        return ssum(lambda n: sg(m + n) / ((2 * m + 1) * (2 * n + 1))
                    * sech(pi * sqrt(((2 * m + 1) / a) ** 2 + ((2 * n + 1) / b) ** 2) * L / 2))
    return c * ssum(row)


def box_sides(a, b, L):
    """1 - P(sides): z-modes cos((2k+1) pi z/L); v_k = centre value of (Lap - mu^2) v = 0 on the
    a x b rectangle with v = 1 on its boundary (two wall pairs, each a 1-D series)."""
    def v(k):
        mu = (2 * k + 1) * pi / L
        wx = 4 / pi * ssum(lambda n: sg(n) / (2 * n + 1) * sech(sqrt(((2 * n + 1) * pi / b) ** 2 + mu ** 2) * a / 2))
        wy = 4 / pi * ssum(lambda m: sg(m) / (2 * m + 1) * sech(sqrt(((2 * m + 1) * pi / a) ** 2 + mu ** 2) * b / 2))
        return wx + wy
    return 4 / pi * ssum(lambda k: sg(k) / (2 * k + 1) * v(k), tol=mpf(10) ** -122)
show('bc10.4', box_ends(1, 2, 3), 1 - box_sides(1, 2, 3))

# bc10.5  semi-infinite strip {x > 0, 0 < y < 1}, start (1, 1/3), P(end x = 0)
A5 = 4 / pi * NS(lambda k: sin((2 * k + 1) * pi / 3) / (2 * k + 1) * exp(-(2 * k + 1) * pi))
show('bc10.5', A5, 2 / pi * atan(sin(pi / 3) / sinh(pi)))

# bc10.6  10 x 1 rectangle, drift: generator (1/2) Lap + d/dx, centre, P(ends)
c = mpf(1)
A6 = 4 / pi * NS(lambda k: sg(k) / (2 * k + 1) * cosh(5 * c) / cosh(5 * sqrt(c ** 2 + ((2 * k + 1) * pi) ** 2)))
def b6(n):
    n = 2 * n + 1                           # even n vanish (sin(n pi/2))
    om = n * pi / 10
    bn = om / 5 * exp(-5 * c) * (1 + exp(10 * c)) / (c ** 2 + om ** 2)    # (-1)^n = -1
    return bn * sg((n - 1) // 2) * sech(sqrt(om ** 2 + c ** 2) / 2)
show('bc10.6', A6, 1 - ssum(b6))

# bc10.7  sector of angle pi/3, radius 1, start r = 1/2 on the bisector: P(arc)
A7 = 4 / pi * NS(lambda k: sg(k) / (2 * k + 1) * mpf(8) ** -(2 * k + 1))
show('bc10.7', A7, 2 / pi * atan(mpf(16) / 63))

# bc10.8  cylinder radius 1/2, length 10, centre: P(ends)
A8 = nsum(lambda k: 2 / (besseljzero(0, int(k)) * besselj(1, besseljzero(0, int(k))))
          * sech(10 * besseljzero(0, int(k))), [1, inf])
B8 = 1 - 4 / pi * ssum(lambda k: sg(k) / (2 * k + 1) / besseli(0, (2 * k + 1) * pi / 20))
show('bc10.8', A8, B8)

# bc10.9  unit square, start (1/2, 99/100), P(top side)
A9 = 4 / pi * NS(lambda k: sg(k) / (2 * k + 1) * sinh((2 * k + 1) * pi * mpf(99) / 100) / sinh((2 * k + 1) * pi))
B9 = 1 - sq_side(1, 1, mpf(1) / 2, mpf(1) / 100) - 2 * sq_side(1, 1, mpf(99) / 100, mpf(1) / 2)
show('bc10.9', A9, B9)

# bc10.10 box 1 x 1 x 10, centre, P(ends)
show('bc10.10', box_ends(1, 1, 10), 1 - box_sides(1, 1, 10))
