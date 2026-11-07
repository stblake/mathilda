"""Holdout A, family C10 (experiment 205): references for hc10.1 .. hc10.10.

Brownian motion in a box: hitting probabilities u (harmonic, u = 1 on the
target part of the boundary, 0 elsewhere) and an expected exit time w
(Delta w = -2, w = 0 on the boundary, standard BM with generator Delta/2).

Method A: the separated series in the variable the statement uses (the case body).
Method B: an INDEPENDENT identity -- the separated series of the complementary
targets in the OTHER variable(s), using sum of all hitting probabilities = 1
(or, for w, the expansion in the other variable), or a closed form.
Both at dps 90; the script prints both and their relative difference.
Usage: python3.11 hc10_ref.py
"""
from mpmath import mp, mpf, pi, cos, sin, exp, sech, nsum, inf, nstr, atan, sinh, sqrt, fsum

mp.dps = 90


def ratio_cosh(c, x, a):
    """cosh(c x)/cosh(c a), stable."""
    return (exp(-c*(a-x)) + exp(-c*(a+x)))/(1 + exp(-2*c*a))


def uE(a, b, x, y):
    """rect [-a,a]x[-b,b], u = 1 on x = +-a."""
    return 4/pi*nsum(lambda k: (-1)**int(k)/(2*k+1)*cos((2*k+1)*pi*y/(2*b))
                     * ratio_cosh((2*k+1)*pi/(2*b), x, a), [0, inf])


def uS(a, b, x, y):
    """u = 1 on y = +-b: the series in x.  Its terms decay only like
    exp(-(2k+1) pi (b-|y|)/(2a)) with a non-alternating cos(.. x) pattern, which
    defeats nsum's extrapolation, so it is summed directly (4000 terms: the
    tail is < 1e-100 for every case here)."""
    return 4/pi*fsum((-1)**k/mpf(2*k+1)*cos((2*k+1)*pi*x/(2*a))
                     * ratio_cosh((2*k+1)*pi/(2*a), y, b) for k in range(4000))


def sq_bottom(x, y):
    """unit square [0,1]^2, u = 1 on y = 0."""
    return 4/pi*nsum(lambda j: sin((2*j+1)*pi*x)/(2*j+1)
                     * (exp(-(2*j+1)*pi*y) - exp(-(2*j+1)*pi*(2-y)))/(1 - exp(-2*(2*j+1)*pi)), [0, inf])


def wE(a, b, x, y):
    """expected exit time, rect [-a,a]x[-b,b]: (b^2-y^2) - sum A_k cos(c y) cosh(c x)/cosh(c a)."""
    return (b*b - y*y) - nsum(lambda k: 32*b*b*(-1)**int(k)/((2*k+1)**3*pi**3)
                              * cos((2*k+1)*pi*y/(2*b)) * ratio_cosh((2*k+1)*pi/(2*b), x, a), [0, inf])


def box_ends(a, b, c):
    """box [-a,a]x[-b,b]x[-c,c], centre, u = 1 on x = +-a."""
    return 16/pi**2*nsum(lambda j, k: (-1)**int(j+k)/((2*j+1)*(2*k+1))
                         * sech(a*pi*sqrt(((2*j+1)/(2*b))**2 + ((2*k+1)/(2*c))**2)), [0, inf], [0, inf])


def cube_face(x, y, z):
    """unit cube, u = 1 on z = 0: sum_{m,n odd} 16/(pi^2 m n) sin(m pi x) sin(n pi y)
    sinh(pi r (1-z))/sinh(pi r), r = sqrt(m^2+n^2)."""
    def t(j, k):
        m, n = 2*j+1, 2*k+1
        r = pi*sqrt(m*m + n*n)
        return 16/(pi**2*m*n)*sin(m*pi*x)*sin(n*pi*y)*(exp(-r*z) - exp(-r*(2-z)))/(1 - exp(-2*r))
    return nsum(t, [0, inf], [0, inf])


h = mpf(1)/2
cases = []
# 1: 4x1 rectangle, centre, ends
cases.append(("hc10.1", uE(2, h, 0, 0), 1 - uS(2, h, 0, 0)))
# 2: 10x1, start (3,0), ends
cases.append(("hc10.2", uE(5, h, 3, 0), 1 - uS(5, h, 3, 0)))
# 3: 10x1, start (4, 1/4), ends
cases.append(("hc10.3", uE(5, h, 4, mpf(1)/4), 1 - uS(5, h, 4, mpf(1)/4)))
# 4: 2x1 rectangle, centre, the LONG sides y = +-1/2
cases.append(("hc10.4", uS(1, h, 0, 0), 1 - uE(1, h, 0, 0)))
# 5: unit square, start (1/2, 1/4), bottom side only
x5, y5 = h, mpf(1)/4
B5 = 1 - (sq_bottom(x5, 1-y5) + sq_bottom(y5, x5) + sq_bottom(y5, 1-x5))   # top, left, right by rotation
cases.append(("hc10.5", sq_bottom(x5, y5), B5))
# 6: semi-infinite strip [0,inf)x[0,1], start (1,1/2), hit the end x = 0
cases.append(("hc10.6", 4/pi*nsum(lambda k: (-1)**int(k)*exp(-(2*k+1)*pi)/(2*k+1), [0, inf]),
              2/pi*atan(1/sinh(pi))))
# 7: expected exit time, 2x1 rectangle [-1,1]x[-1/2,1/2], centre
cases.append(("hc10.7", wE(1, h, 0, 0), wE(h, 1, 0, 0)))
# 8: box 3x1x1, centre, the two end faces
B8 = 1 - 2*box_ends(h, mpf(3)/2, h)
cases.append(("hc10.8", box_ends(mpf(3)/2, h, h), B8))
# 9: unit cube, start (1/2,1/2,1/4), bottom face
# complement: top face (z -> 1-z) + four side faces; side x=0 at (1/2,1/2,1/4) is
# cube_face with the roles (x,y,z) -> (y, z, x) = (1/2, 1/4, 1/2); all four sides equal by symmetry
B9 = 1 - cube_face(h, h, 1 - mpf(1)/4) - 4*cube_face(h, mpf(1)/4, h)
cases.append(("hc10.9", cube_face(h, h, mpf(1)/4), B9))
# 10: 4x1 rectangle, start (1.9, 0) near an end
x10 = mpf(19)/10
cases.append(("hc10.10", uE(2, h, x10, 0), 1 - uS(2, h, x10, 0)))

for lab, A, B in cases:
    print(lab, nstr(A, 70), nstr(B, 70), "reldiff", nstr(abs(A-B)/abs(A), 3))
