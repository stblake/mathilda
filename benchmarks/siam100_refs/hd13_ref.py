"""Holdout A, family d13 (experiment 218): method 1 references (closed forms,
mpmath dps 80; hd13.6 and hd13.9 by mpmath.odefun shooting).  Method 2 is
hd13_ref.wl (NDSolve WorkingPrecision 80-90 + FindRoot).

  1  y' = y(1-y), y(10) = 1/2:                  y0 = 1/(1+e^10)
  2  u' = A u, A = [[-1,10,0],[0,-1,10],[0,0,-1]], u(0) = (0,0,z0), |u(5)| = 1:
                                                z0 = e^5/Sqrt[2500*5^4 + 100*5^2 + 1]
  3  y' = y^2 - y, y(3) = 2:                    y0 = 1/(1 - e^-3/2)
  4  th'' = -sin th, th(0)=0, th'(0)=v0, th(2) = 1:  th = 2 asin(k sn(t|k^2)), k = v0/2
  5  y'' = -e^y, y(0)=0, y'(0)=s, y(1) = 0, smaller s:  s = th tanh(th/4), th = Sqrt2 cosh(th/4) (smaller root)
  6  f''' + f f''/2 = 0, f(0)=f'(0)=0, f''(0)=s, f'(10) = 1      (odefun shooting)
  7  r'' = -r/|r|^3, r(0)=(1,0), r'(0)=(0,v0), |r(2)| = 3/2:  Kepler's equation
  8  y' = t^2 + y^2, y(0)=y0, y(1) = 2:        u'' + t^2 u = 0, y = -u'/u, u = 0F1(;3/4;-t^4/16) - y0 t 0F1(;5/4;-t^4/16)
  9  exp 188's ODE, u(0) = (0,0,z0), |u(30)| = 1  (odefun shooting, secant in the needle window)
 10  x'' + x + x^3 = 0, x(0)=x0, x'(0)=0, x(10) = 0, smallest x0 > 0:
                                                x = x0 cn(w t|m), w^2 = 1+x0^2, m = x0^2/(2(1+x0^2)); 10 w = 7 K(m)
"""
import sys
import mpmath as mp
mp.mp.dps = 80
R = {}
R[1] = 1/(1 + mp.e**10)
R[2] = mp.e**5/mp.sqrt(2500*mp.mpf(5)**4 + 100*25 + 1)
R[3] = 1/(1 - mp.e**-3/2)
R[4] = 2*mp.findroot(lambda k: k*mp.ellipfun('sn', 2, m=k*k) - mp.sin(mp.mpf(1)/2), mp.mpf("0.55"))
th5 = mp.findroot(lambda th: th - mp.sqrt(2)*mp.cosh(th/4), mp.mpf("1.5"))
R[5] = th5*mp.tanh(th5/4)
def kepler_r(v0, t):
    a = 1/(2 - v0*v0); e = v0*v0 - 1; n = a**mp.mpf(-1.5)
    M = n*t
    E = mp.findroot(lambda E: E - e*mp.sin(E) - M, M)
    return a*(1 - e*mp.cos(E))
R[7] = mp.re(mp.findroot(lambda v: kepler_r(v, 2) - mp.mpf(3)/2, mp.mpf("1.2")))
u1 = lambda t: mp.hyp0f1(mp.mpf(3)/4, -t**4/16)
u2 = lambda t: t*mp.hyp0f1(mp.mpf(5)/4, -t**4/16)
one = mp.mpf(1)
R[8] = (2*u1(one) + mp.diff(u1, one))/(mp.diff(u2, one) + 2*u2(one))
def km(x0):
    w = mp.sqrt(1 + x0*x0); m = x0*x0/(2*(1 + x0*x0))
    return 10*w - 7*mp.ellipk(m)
R[10] = mp.findroot(km, mp.mpf("0.45"))
if "fast" not in sys.argv:
    def blas(s):
        return mp.odefun(lambda t, u: [u[1], u[2], -u[0]*u[2]/2], 0, [mp.mpf(0), mp.mpf(0), s])(10)[1] - 1
    R[6] = mp.findroot(blas, (mp.mpf("0.332057"), mp.mpf("0.332058")))
    a7, a03, a1 = mp.mpf(7)/10, mp.mpf(3)/100, mp.mpf(1)/10
    def F9(t, u):
        x, y, z = u
        n = mp.sqrt(x*x + y*y + z*z)
        return [-x + 10*y + n*(-a7*y - a03*z), -y + 10*z + n*(a7*y - a1*z), -z + n*(a03*x + a1*y)]
    def n30(z0):
        u = mp.odefun(F9, 0, [mp.mpf(0), mp.mpf(0), z0])(30)
        return mp.sqrt(sum(c*c for c in u)) - 1
    R[9] = mp.findroot(n30, (mp.mpf("0.0064839843787045"), mp.mpf("0.0064839843787046")))
for k in sorted(R):
    print("hd13.%d" % k, mp.nstr(R[k], 70))
