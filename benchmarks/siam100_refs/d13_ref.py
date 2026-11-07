"""Reference for SIAM100 Appendix D problem 13 (experiment 188), method A:
mpmath.odefun (adaptive Taylor series) at dps 85, secant (mpmath.findroot) on
z0 for ||u(50)|| = 1.  Method B: d13_ref.wl (NDSolve "Extrapolation", WP 90).

Landscape (scipy DOP853 rtol 1e-13 scan, z0 in [1e-7, 100]): u = 0 is stable;
there is a saddle equilibrium with ||u*|| = 0.21523 (one unstable eigenvalue
0.2995) whose stable manifold meets the z-axis at z_c = 0.0064823961796...;
below z_c the orbit decays (max ||u|| < 0.2153 < 1), above it the orbit leaves
the saddle and blows up in finite time (before t = 50 once z0 > ~z_c(1+1e-9)).
So ||u(50)|| = 1 has exactly ONE solution, on the blow-up side, a hair above z_c.
"""
import mpmath as mp
import time
mp.mp.dps = 85
a7, a03, a1 = mp.mpf(7)/10, mp.mpf(3)/100, mp.mpf(1)/10
def F(t, u):
    x, y, z = u
    n = mp.sqrt(x*x + y*y + z*z)
    return [-x + 10*y + n*(-a7*y - a03*z), -y + 10*z + n*(a7*y - a1*z), -z + n*(a03*x + a1*y)]
def g(z0):
    t0 = time.time()
    sol = mp.odefun(F, 0, [mp.mpf(0), mp.mpf(0), z0])
    u = sol(50)
    r = mp.sqrt(sum(c*c for c in u)) - 1
    print("  eval", mp.nstr(z0, 25), mp.nstr(r, 10), "%.1fs" % (time.time()-t0), flush=True)
    return r
z = mp.findroot(g, (mp.mpf('0.006482396179622983'), mp.mpf('0.006482396179622990')), tol=mp.mpf(10)**-150)
print("z0 =", mp.nstr(z, 75))
