"""Holdout A, family d11 (experiment 216): method 1 references, mpmath dps 80
(closed forms / Kepler's equation / energy quadrature; hd11.8 and hd11.10 by
mpmath.odefun shooting in (theta, T)).  Method 2: hd11_ref.wl (NDSolve
"Extrapolation" WorkingPrecision 80 + FindRoot).  Shortest branches: closed
forms (1-4, 9), hd11_scan.py (5-8, 10).

  1  r'' = (0,-1), r(0) = 0, |r'(0)| = 2, target (1, 1/2):       T^4 - 14T^2 + 5 = 0, T = Sqrt[7 - 2 Sqrt[11]]
  2  r'' = -r, r(0) = (1,0), |r'(0)| = 1, target (0, 1/2):       cos 2T = -1/4
  3  x'' = -x, y'' = -4y, r(0) = 0, |r'(0)| = 1, target (1/2, 1/4):  1/(4 sin^2 T) + 1/(4 sin^2 2T) = 1
  4  r'' = -r'/2 + (0,-1), r(0) = 0, |r'(0)| = 2, target (1, 0): (1/(4E))^2 + (T/(2E) - 1)^2 = 1, E = 1 - e^{-T/2}
  5  r'' = -r/|r|^3, r(0) = (1,0), |r'(0)| = 1, target (-1/2, 1/2)  (Kepler, f and g functions)
  6  r'' = -r/|r|^3, r(0) = (1,0), |r'(0)| = 6/5, target (0, 3/2)   (Kepler)
  7  unit masses at (2, +-1), r(0) = 0, |r'(0)| = 1, target (4, 0): theta = 0 by symmetry,
     T = Int_0^4 dx/Sqrt[1 + 2(U(x) - U(0))], U(x) = 2/Sqrt[(x-2)^2 + 1]
  8  unit masses at (1,1), (2,-1), (3,1), r(0) = 0, |r'(0)| = 1, target (4, 0)   (odefun shooting)
  9  r'' = (y', -x') (unit magnetic field), r(0) = 0, |r'(0)| = 1, target (1, 1/2):  T = 2 ArcSin[Sqrt[5]/4]
 10  the six masses of exp 186, r(0) = 0, |r'(0)| = 1, target (1, 2)   (odefun shooting)
"""
import sys
import mpmath as mp
mp.mp.dps = 80
R = {}
R[1] = mp.sqrt(7 - 2*mp.sqrt(11))
R[2] = mp.acos(-mp.mpf(1)/4)/2
R[3] = mp.findroot(lambda t: 1/(4*mp.sin(t)**2) + 1/(4*mp.sin(2*t)**2) - 1, mp.mpf("0.6"))
def f4(t):
    E = 1 - mp.e**(-t/2)
    return (1/(4*E))**2 + (t/(2*E) - 1)**2 - 1
R[4] = mp.findroot(f4, mp.mpf("0.6"))
def kepler(r0, v0, dt):
    """position at time dt on the Kepler ellipse through (r0, v0), mu = 1 (f and g functions)."""
    r = mp.sqrt(r0[0]**2 + r0[1]**2); v2 = v0[0]**2 + v0[1]**2
    a = 1/(2/r - v2); n = a**mp.mpf(-1.5); s0 = r0[0]*v0[0] + r0[1]*v0[1]
    dE = mp.findroot(lambda E: E - (1 - r/a)*mp.sin(E) + s0/mp.sqrt(a)*(1 - mp.cos(E)) - n*dt, n*dt)
    f = 1 - a/r*(1 - mp.cos(dE)); g = dt + (mp.sin(dE) - dE)/n
    return [f*r0[0] + g*v0[0], f*r0[1] + g*v0[1]]
def kep_case(speed, tgt, guess):
    F = lambda th, T: [c - d for c, d in zip(kepler([mp.mpf(1), mp.mpf(0)], [speed*mp.cos(th), speed*mp.sin(th)], T), tgt)]
    return mp.findroot(F, guess)
S5 = kep_case(mp.mpf(1), [-mp.mpf(1)/2, mp.mpf(1)/2], (mp.mpf("1.887059167646"), mp.mpf("1.498227157579")))
R[5] = S5[1]
S6 = kep_case(mp.mpf(6)/5, [mp.mpf(0), mp.mpf(3)/2], (mp.mpf("1.542468565154"), mp.mpf("1.804739896993")))
R[6] = S6[1]
U = lambda x: 2/mp.sqrt((x - 2)**2 + 1)
R[7] = mp.quad(lambda x: 1/mp.sqrt(1 + 2*(U(x) - U(0))), [0, 2, 4])
R[9] = 2*mp.asin(mp.sqrt(5)/4)
def masses_case(ms, tgt, guess):
    def rhs(t, u):
        x, y, vx, vy = u
        ax = ay = mp.mpf(0)
        for (a, b) in ms:
            dx, dy = a - x, b - y
            r3 = (dx*dx + dy*dy)**mp.mpf(1.5)
            ax += dx/r3; ay += dy/r3
        return [vx, vy, ax, ay]
    cache = {}
    def F(th, T):
        if th not in cache:
            cache[th] = mp.odefun(rhs, 0, [mp.mpf(0), mp.mpf(0), mp.cos(th), mp.sin(th)])
        u = cache[th](T)
        return [u[0] - tgt[0], u[1] - tgt[1]]
    return mp.findroot(F, guess, tol=mp.mpf(10)**-140)
if "fast" not in sys.argv:
    R[8] = masses_case([(1, 1), (2, -1), (3, 1)], (4, 0), (mp.mpf("-0.350867468163"), mp.mpf("2.878059641701")))[1]
    R[10] = masses_case([(2, -1), (2, 0), (2, 1), (3, -1), (3, 0), (3, 1)], (1, 2), (mp.mpf("0.944671460204"), mp.mpf("2.028352767943")))[1]
# launch angles of the closed-form cases (for the starting guesses)
th3 = mp.acos(1/(2*mp.sin(R[3])))
E4 = 1 - mp.e**(-R[4]/2); th4 = mp.atan2(R[4]/(2*E4) - 1, 1/(4*E4))
# magnetic: r(t) = (sin th - sin(th - t), cos(th - t) - cos th); solve r(T) = (1, 1/2) for th
th9 = mp.findroot(lambda th: mp.sin(th) - mp.sin(th - R[9]) - 1, mp.mpf("1"))
print("theta3", mp.nstr(th3, 20), "theta4", mp.nstr(th4, 20), "theta9", mp.nstr(th9, 20),
      "check9", mp.nstr(mp.cos(th9 - R[9]) - mp.cos(th9) - mp.mpf(1)/2, 5))
print("theta5", mp.nstr(S5[0], 30), "theta6", mp.nstr(S6[0], 30))
for k in sorted(R):
    print("hd11.%d" % k, mp.nstr(R[k], 70))
