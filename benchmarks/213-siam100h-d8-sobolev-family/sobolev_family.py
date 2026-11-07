"""Experiment 213 -- SIAM 100-digit HOLDOUT A, family D8 (Python column).

The ten integrals over R^N of sobolev_family.m, as written: scipy.integrate
quad / dblquad / tplquad (epsabs 0, epsrel 1e-13) at machine precision;
mpmath.quad over [-inf, inf]^N (tanh-sinh, N <= 3) at 60 digits."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import vrun, require
import math
import mpmath as M
from mpmath import mpf
from scipy.integrate import quad, dblquad, tplquad

require(["scipy.integrate:quad", "scipy.integrate:dblquad", "scipy.integrate:tplquad", "mpmath:quad"])
inf, pi = math.inf, math.pi
def mach(f, N):
    if N == 1:
        I = quad(f, -inf, inf, epsabs=0, epsrel=1e-13, limit=200)[0]
    elif N == 2:
        I = dblquad(lambda y, x: f(x, y), -inf, inf, -inf, inf, epsabs=0, epsrel=1e-13)[0]
    else:
        I = tplquad(lambda z, y, x: f(x, y, z), -inf, inf, -inf, inf, -inf, inf, epsabs=0, epsrel=1e-13)[0]
    return math.sqrt(I / (2 * pi) ** N)
def wp(f, N):
    return M.sqrt(M.quad(f, *([[-M.inf, M.inf]] * N)) / (2 * M.pi) ** N)

# hd8.1  c_1(R^1)
vrun("hd8.1 machine", lambda: mach(lambda x: 1 / (1 + x**2), 1), 10, 3, cap=30)
vrun("hd8.1 wp50", lambda: wp(lambda x: 1 / (1 + x**2), 1), 48, cap=60, dps=60)
# hd8.2  c_2(R^1)
vrun("hd8.2 machine", lambda: mach(lambda x: 1 / (1 + x**2 + x**4), 1), 10, 3, cap=30)
vrun("hd8.2 wp50", lambda: wp(lambda x: 1 / (1 + x**2 + x**4), 1), 48, cap=60, dps=60)
# hd8.3  c_5(R^1)
vrun("hd8.3 machine", lambda: mach(lambda x: 1 / (1 + x**2 + x**4 + x**6 + x**8 + x**10), 1), 10, 3, cap=30)
vrun("hd8.3 wp50", lambda: wp(lambda x: 1 / (1 + x**2 + x**4 + x**6 + x**8 + x**10), 1), 48, cap=60, dps=60)
# hd8.4  c_10(R^1)
vrun("hd8.4 machine", lambda: mach(lambda x: 1 / (1 + x**2 + x**4 + x**6 + x**8 + x**10 + x**12 + x**14 + x**16 + x**18 + x**20), 1), 10, 3, cap=30)
vrun("hd8.4 wp50", lambda: wp(lambda x: 1 / (1 + x**2 + x**4 + x**6 + x**8 + x**10 + x**12 + x**14 + x**16 + x**18 + x**20), 1), 48, cap=60, dps=60)
# hd8.5  c_3(R^2)
vrun("hd8.5 machine", lambda: mach(lambda x, y: 1 / (1 + x**2 + y**2 + x**4 + x**2*y**2 + y**4 + x**6 + x**4*y**2 + x**2*y**4 + y**6), 2), 10, 3, cap=30)
vrun("hd8.5 wp50", lambda: wp(lambda x, y: 1 / (1 + x**2 + y**2 + x**4 + x**2*y**2 + y**4 + x**6 + x**4*y**2 + x**2*y**4 + y**6), 2), 48, cap=60, dps=60)
# hd8.6  c_4(R^2)
vrun("hd8.6 machine", lambda: mach(lambda x, y: 1 / (1 + x**2 + y**2 + x**4 + x**2*y**2 + y**4 + x**6 + x**4*y**2 + x**2*y**4 + y**6 + x**8 + x**6*y**2 + x**4*y**4 + x**2*y**6 + y**8), 2), 10, 3, cap=30)
vrun("hd8.6 wp50", lambda: wp(lambda x, y: 1 / (1 + x**2 + y**2 + x**4 + x**2*y**2 + y**4 + x**6 + x**4*y**2 + x**2*y**4 + y**6 + x**8 + x**6*y**2 + x**4*y**4 + x**2*y**6 + y**8), 2), 48, cap=60, dps=60)
# hd8.7  c_2(R^3)
vrun("hd8.7 machine", lambda: mach(lambda x, y, z: 1 / (1 + x**2 + y**2 + z**2 + x**4 + x**2*y**2 + x**2*z**2 + y**4 + y**2*z**2 + z**4), 3), 10, 3, cap=30)
vrun("hd8.7 wp50", lambda: wp(lambda x, y, z: 1 / (1 + x**2 + y**2 + z**2 + x**4 + x**2*y**2 + x**2*z**2 + y**4 + y**2*z**2 + z**4), 3), 48, cap=60, dps=60)
# hd8.8  Bessel potential (1+|x|^2)^(-3/2) on R^2
vrun("hd8.8 machine", lambda: mach(lambda x, y: 1 / ((1 + x*x + y*y)**1.5), 2), 10, 3, cap=30)
vrun("hd8.8 wp50", lambda: wp(lambda x, y: 1 / ((1 + x*x + y*y)**mpf(1.5)), 2), 48, cap=60, dps=60)
# hd8.9  Bessel potential (1+|x|^2)^(-5/2) on R^3
vrun("hd8.9 machine", lambda: mach(lambda x, y, z: 1 / ((1 + x*x + y*y + z*z)**2.5), 3), 10, 3, cap=30)
vrun("hd8.9 wp50", lambda: wp(lambda x, y, z: 1 / ((1 + x*x + y*y + z*z)**mpf(2.5)), 3), 48, cap=60, dps=60)
# hd8.10  anisotropic 1/(1+x^2+y^4) on R^2
vrun("hd8.10 machine", lambda: mach(lambda x, y: 1 / (1 + x*x + y**4), 2), 10, 3, cap=30)
vrun("hd8.10 wp50", lambda: wp(lambda x, y: 1 / (1 + x*x + y**4), 2), 48, cap=60, dps=60)
