"""Holdout A, family d15 (experiment 220): method 1 references, mpmath dps 80.
Method 2 for each is printed alongside (a second formula / quadrature / ODE
solve, mpmath) or is hd15_ref.wl (BesselJZero, NIntegrate).

  1  y' = 1 + y^2, y(0) = 0:                    T = pi/2
  2  y' = t^2 + y^2, y(0) = 0:                  T = Sqrt[2 j], j the first zero of J_{-1/4}
       (A: findroot on besselj(-1/4, t^2/2); B: findroot on 0F1(;3/4;-t^4/16); C: WL BesselJZero)
  3  y'' = y^2, y(0)=1, y'(0)=0:                T = Sqrt[3/2] B(1/6,1/2)/3   (B: quad of the energy integral)
  4  y' = e^y - y, y(0) = 0:                    T = Int_0^inf dy/(e^y - y)  (A: tanh-sinh; B: Gauss-Legendre; C: WL NIntegrate)
  5  x' = x^2 + y^2, y' = 2xy, (1, 1/2):        u = x+y has u' = u^2, u0 = 3/2: T = 2/3 (x-y blows up only at t = 2)
  6  y' = y^2 + e^t, y(0) = 0:                  y = -u'/u, u'' + e^t u = 0, u(0)=1, u'(0)=0; first zero of u
       (A: u = a J0(2e^{t/2}) + b Y0(2e^{t/2}); B: odefun on the linear ODE for u)
  7  u_t = u_xx + u + u^2 on [0,1], u_x = 0 at both ends, u(x,0) = 1:  spatially constant, T = log 2
  8  u_t = Lap u + e^u on [0,3]^2, du/dn = 0, u(x,y,0) = 0:          spatially constant, T = 1
  9  u_t = u_xx, u(0,t) = u(1,t) = 0, u(x,0) = 4x(1-x); time at which u(1/2,t) = 1/2
       (A: sine series; B: method of images, Gaussian heat kernel, mpmath quad)
 10  y'' = e^y, y(0)=y'(0)=0:                   T = pi/Sqrt[2]   (B: quad of Int_0^inf dy/Sqrt[2(e^y-1)])
"""
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi
A, B = {}, {}
A[1] = pi/2; B[1] = mp.quad(lambda y: 1/(1 + y*y), [0, mp.inf])
A[2] = mp.findroot(lambda t: mp.besselj(-mp.mpf(1)/4, t*t/2), mp.mpf(2))
B[2] = mp.findroot(lambda t: mp.hyp0f1(mp.mpf(3)/4, -t**4/16), mp.mpf(2))
A[3] = mp.sqrt(mp.mpf(3)/2)*mp.beta(mp.mpf(1)/6, mp.mpf(1)/2)/3
B[3] = mp.quad(lambda w: 2/mp.sqrt(2*(1 - w**6)/3), [0, mp.mpf(1)/2]) + mp.quad(lambda v: 4/mp.sqrt(2*(lambda s: 6 - 15*s + 20*s**2 - 15*s**3 + 6*s**4 - s**5)(v*v)/3), [0, mp.sqrt(mp.mpf(1)/2)])   # y = w^-2, w = 1 - v^2
A[4] = mp.quad(lambda y: 1/(mp.e**y - y), [0, 1, 10, mp.inf])
B[4] = mp.quad(lambda y: 1/(mp.e**y - y), [0, 1, 10, 60, 200], method="gauss-legendre") \
     + mp.quad(lambda y: 1/(mp.e**y - y), [200, mp.inf])
A[5] = mp.mpf(2)/3; B[5] = mp.quad(lambda u: 1/u**2, [mp.mpf(3)/2, mp.inf])
J0, Y0, J1, Y1 = [f(2) for f in (lambda s: mp.besselj(0, s), lambda s: mp.bessely(0, s),
                                  lambda s: mp.besselj(1, s), lambda s: mp.bessely(1, s))]
# a J0(2) + b Y0(2) = 1,  a J1(2) + b Y1(2) = 0
det = J0*Y1 - Y0*J1
a6, b6 = Y1/det, -J1/det
u6 = lambda t: a6*mp.besselj(0, 2*mp.e**(t/2)) + b6*mp.bessely(0, 2*mp.e**(t/2))
A[6] = mp.findroot(u6, mp.mpf("1.0"))
sol6 = mp.odefun(lambda t, u: [u[1], -mp.e**t*u[0]], 0, [mp.mpf(1), mp.mpf(0)])
B[6] = mp.findroot(lambda t: sol6(t)[0], mp.mpf("1.0"))
A[7] = mp.log(2); B[7] = mp.quad(lambda u: 1/(u + u*u), [1, mp.inf])
A[8] = mp.mpf(1); B[8] = mp.quad(lambda u: mp.e**-u, [0, mp.inf])
def ser9(t):
    s = 0
    for n in range(1, 400, 2):
        term = 32/(n**3*pi**3)*(-1)**((n - 1)//2)*mp.e**(-n*n*pi*pi*t)
        s += term
        if abs(term) < mp.mpf(10)**-90: break
    return s
def img9(t):
    g = lambda z: mp.e**(-z*z/(4*t))/mp.sqrt(4*pi*t)
    x = mp.mpf(1)/2
    return mp.quad(lambda xi: 4*xi*(1 - xi)*sum(g(x - xi + 2*k) - g(x + xi + 2*k) for k in range(-6, 7)), [0, x, 1])
A[9] = mp.findroot(lambda t: ser9(t) - mp.mpf(1)/2, mp.mpf("0.07"))
B[9] = mp.findroot(lambda t: img9(t) - mp.mpf(1)/2, A[9]*(1 + mp.mpf(10)**-20))
A[10] = pi/mp.sqrt(2); B[10] = mp.quad(lambda s: 2*s/mp.sqrt(2*mp.expm1(s*s)), [0, 1]) + mp.quad(lambda y: 1/mp.sqrt(2*mp.expm1(y)), [1, 10, 100, 300]) + mp.quad(lambda y: 1/mp.sqrt(2*mp.expm1(y)), [300, mp.inf])
for k in range(1, 11):
    print("hd15.%d" % k, mp.nstr(A[k], 70), " agree:", int(-mp.log10(abs(A[k] - B[k])/A[k] + mp.mpf(10)**-80)))
