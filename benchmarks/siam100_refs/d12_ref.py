"""Reference for SIAM100 Appendix D problem 12 (experiment 187).

Limiting (eps -> 0) dynamics (Rubin-Ungar / Takens / Bornemann, "Homogenization
in time of singularly perturbed mechanical systems", LNM 1687): the particle
moves on the unit circle, theta = polar angle, under the potential
    V(theta) = sin(theta) + theta* omega(x),   x = cos(theta),
where omega(x) = 2 sqrt(1 + alpha x^2) is the normal frequency of the constraint
(Hessian of (1/2)(1+alpha x^2)(r^2-1)^2 in r at r = 1), and theta* = E_perp(0)/omega(x0)
= (1/2)/2 = 1/4 is the adiabatic invariant of the normal oscillation (initial
normal velocity 1, x0 = 0).  So V = sin(theta) + (1/2) sqrt(1 + alpha cos^2 theta),
theta(0) = pi/2, theta'(0) = -1, energy E = 1/2 + V(pi/2) = 2.
(Checked against the full eps-system with eps = 1e-2 .. 2.5e-3: hitting times
agree to O(eps^2); the naive constraint without the correction is off by 3-15%.)

Method A (this file): first hitting time of y = 0 by energy quadrature
    T(alpha) = Int_0^{pi/2} dtheta / sqrt(2 (E - V(theta)))
and secant root of T(alpha) = 10, mpmath dps 80.
Method B: d12_ref.wl (NDSolve WorkingPrecision 80 + WhenEvent + FindRoot).
"""
import mpmath as mp
mp.mp.dps = 80
def V(th, a): return mp.sin(th) + mp.sqrt(1 + a*mp.cos(th)**2)/2
def thmax(a):   # barrier: alpha sin = 2 sqrt(1 + alpha cos^2)  (alpha > 2)
    return mp.findroot(lambda th: a*mp.sin(th) - 2*mp.sqrt(1 + a*mp.cos(th)**2), mp.mpf(1))
def T(a):
    tm = thmax(a)
    f = lambda th: 1/mp.sqrt(2*(2 - V(th, a)))
    return mp.quad(f, mp.linspace(0, tm, 6) + mp.linspace(tm, mp.pi/2, 6)[1:])
# critical alpha: max V = 2
def crit(a):
    return V(thmax(a), a) - 2
ac = mp.findroot(crit, mp.mpf(8))
print("alpha_c =", mp.nstr(ac, 30))
a = mp.findroot(lambda a: T(a) - 10, (ac - mp.mpf('1e-6'), ac - mp.mpf('2e-6')))
print("alpha =", mp.nstr(a, 70))
print("T check", mp.nstr(T(a), 30))
