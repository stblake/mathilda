"""Holdout A, family d12 (experiment 217): method 1 references (closed forms /
energy quadratures, mpmath dps 80).  Method 2 is hd12_ref.wl (NDSolve +
WhenEvent + FindRoot at WorkingPrecision 80); for hd12.8 and hd12.10, which
have no closed form, method 1 is mpmath.odefun (Taylor) fixed-time shooting.

  1  y'' = -1 - a y'^2, y(0)=0, y'(0)=2; first y'=0 at t=1:  atan(2 sqrt a)/sqrt a = 1
  2  th'' = -a sin th, th(0)=0, th'(0)=1; first th'=0 at t=2:  K(1/(4a))/sqrt a = 2
  3  x'' = -x - a x^3, x(0)=1, x'(0)=0; first x=0 at t=6/5:  K(a/(2(1+a)))/sqrt(1+a) = 6/5
  4  x'' + a x' + x = 0, x(0)=1, x'(0)=0; first x=0 at t=2:  w t = pi - atan(2w/a), w = sqrt(1-a^2/4)
  5  y' = a + y^2, y(0)=0; first y=1 at t=1/2:  atan(1/sqrt a)/sqrt a = 1/2
  6  r'' = -a/r^2, r(0)=1, r'(0)=0; first r=1/2 at t=1:  a = (1/2 + pi/4)^2/2
  7  y'' + a t y = 0, y(0)=1, y'(0)=0; first y=0 at t=1:  a = s1^3, s1 the first zero of Bi'(0)Ai(-s) - Ai'(0)Bi(-s)
  8  y'' + (a + cos t) y = 0, y(0)=1, y'(0)=0; first y=0 at t=1   (odefun shooting)
  9  th'' = -cos th + a sin th cos th/(2 sqrt(1 + a cos^2 th)), th(0)=pi/2, th'(0)=-1;
     first sin th = 0 at t=5  (energy quadrature, as siam100_refs/d12_ref.py)
 10  x'' - a(1-x^2)x' + x = 0, x(0)=2, x'(0)=0; first x=0 at t=2   (odefun shooting)
"""
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi
R = {}
R[1] = mp.findroot(lambda a: mp.atan(2*mp.sqrt(a))/mp.sqrt(a) - 1, mp.mpf("1.36"))
R[2] = mp.findroot(lambda a: mp.ellipk(1/(4*a))/mp.sqrt(a) - 2, mp.mpf("0.7"))
R[3] = mp.findroot(lambda a: mp.ellipk(a/(2*(1+a)))/mp.sqrt(1+a) - mp.mpf(6)/5, mp.mpf("0.6"))
def f4(a):
    w = mp.sqrt(1 - a*a/4)
    return (pi - mp.atan(2*w/a))/w - 2
R[4] = mp.findroot(f4, mp.mpf("0.5"))
R[5] = mp.findroot(lambda a: mp.atan(1/mp.sqrt(a))/mp.sqrt(a) - mp.mpf(1)/2, mp.mpf("2.5"))
R[6] = (mp.mpf(1)/2 + pi/4)**2/2
g7 = lambda s: mp.airybi(0, 1)*mp.airyai(-s) - mp.airyai(0, 1)*mp.airybi(-s)
R[7] = mp.findroot(g7, mp.mpf("1.9"))**3
def shoot(rhs, u0, T, idx):
    return lambda a: mp.odefun(lambda t, u: rhs(t, u, a), 0, u0)(T)[idx]
R[8] = mp.findroot(shoot(lambda t, u, a: [u[1], -(a + mp.cos(t))*u[0]], [mp.mpf(1), mp.mpf(0)], 1, 0),
                   (mp.mpf("1.4"), mp.mpf("1.5")))
def V(th, a): return mp.sin(th) + mp.sqrt(1 + a*mp.cos(th)**2)/2
def T9(a):
    tm = mp.findroot(lambda th: a*mp.sin(th) - 2*mp.sqrt(1 + a*mp.cos(th)**2), mp.mpf(1))
    f = lambda th: 1/mp.sqrt(2*(2 - V(th, a)))
    return mp.quad(f, mp.linspace(0, tm, 6) + mp.linspace(tm, pi/2, 6)[1:])
import sys
a9 = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else None
R[9] = mp.findroot(lambda a: T9(a) - 5, (mp.mpf("10.57"), mp.mpf("10.575")))
R[10] = mp.findroot(shoot(lambda t, u, a: [u[1], a*(1 - u[0]**2)*u[1] - u[0]], [mp.mpf(2), mp.mpf(0)], 2, 0),
                    (mp.mpf("0.75"), mp.mpf("0.76")))
if __name__ == "__main__":
    for k in sorted(R):
        print("hd12.%d" % k, mp.nstr(R[k], 70))
