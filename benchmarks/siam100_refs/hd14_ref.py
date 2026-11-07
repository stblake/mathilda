"""Holdout A, family d14 (experiment 219): references, mpmath dps 70.
Each problem has a closed form or a rapidly convergent series (method A) and
an independent second evaluation (method B); both are printed with the number
of digits they agree on.

  1  -u'' = e^{a x} on [0,1], u(0)=0, u'(1)=0;  u(1) = e^a/4
       A: u(1)e^-a = 1/a - (1-e^-a)/a^2;   B: odefun shooting (linear: u = U + c V)
  2  -u'' + a u = 1 on [0,1], u(0)=0, u'(1)=0;  u(1) = 1/4
       A: (1 - sech Sqrt[a])/a;            B: odefun shooting
  3  -u'' = l e^u on [0,1], u(0)=u(1)=0 (lower branch);  u(1/2) = 1/2
       A: l = th^2 e^{-1/2}/2, th = 4 ArcCosh[e^{1/4}];  B: hd14_ref.wl, NDSolve BVP WorkingPrecision 80 + FindRoot
       (the odefun shooting here, argument b3, is too slow at dps 70)
  4  -Lap u = e^{a r^2} on the unit disk, u = 0 on the circle;  u(0,0) = 1/2
       A: u(0) = Ein(a)/(4a) via ei;   B: the power series Sum a^k/(k k!)/(4a)
  5  -Lap u = 1 on [0,s] x [0,1], u = 0;  Int u = 1/20
       A: s/12 - (16/pi^5) Sum_odd tanh(n pi s/2)/n^5;  B: s^3/12 - (16 s^4/pi^5) Sum_odd tanh(n pi/(2s))/n^5
  6  -Lap u = 1 on [0,1] x [0,b], u = 0 on x=0, x=1, y=0, du/dn = 0 on y=b;  Int_0^1 u(x,b) dx = 1/20
       A: Sum_odd 8/(m^4 pi^4)(1 - sech(m pi b));  B: y-modes sin((k+1/2) pi y/b), alternating series (nsum)
  7  -Lap u = 1 on the equilateral triangle of side s, u = 0;  Int u = 1/1000
       A: Int u = Sqrt[3] s^4/320 (Saint-Venant);  B: hd14_ref.wl (exact torsion function, symbolic Integrate)
  8  -Lap u = 1 on the annulus rho < r < 1, u = 0 on r = 1, du/dn = 0 on r = rho;
     Int_{r=rho} u ds = 1/5, the root in (1/2, 9/10)
       A: 2 pi rho ((1-rho^2)/4 + (rho^2/2) log rho);  B: odefun shooting of the radial ODE
  9  -Lap u + a u = 1 on [0,1]^2, u = 0;  Int u = 1/50
       A: Sum_odd 8/(m^2 pi^2 k^2)(1 - (2/k) tanh(k/2)), k^2 = m^2 pi^2 + a;
       B: the double eigen-series Sum_odd,odd 64/(pi^4 m^2 n^2 (pi^2(m^2+n^2) + a)) (nested nsum)
 10  -Lap u = e^{a x} on [0,1]^2, u = 0 on x=0, y=0, y=1, du/dn = 0 on x=1;  Int_0^1 u(1,y) dy = e^a/10
       A: y-modes sin(n pi y), closed-form x-ODE;  B: x-modes sin((j+1/2) pi x), closed-form y-ODE (nsum)
"""
import sys
import mpmath as mp
mp.mp.dps = 70
pi = mp.pi
A, B = {}, {}
def agree(x, y): return int(-mp.log10(abs(x - y)/abs(x) + mp.mpf(10)**-75))
def shoot_lin(f, a, ucond):
    """-u'' = f(x,u) linear BVP on [0,1] via two odefun IVPs: returns the solution object for slope c."""
    pass

# 1
f1 = lambda a: 1/a - (1 - mp.e**-a)/a**2 - mp.mpf(1)/4
A[1] = mp.findroot(f1, mp.mpf("1.5"))
def g1(a):   # u'' = -e^{ax}, u(0)=0, u'(0)=c: u'(1) = c - (e^a-1)/a = 0 exactly -> c; then u(1) by odefun
    c = mp.quad(lambda x: mp.e**(a*x), [0, 1])
    u = mp.odefun(lambda x, w: [w[1], -mp.e**(a*x)], 0, [mp.mpf(0), c])(1)
    return u[0]*mp.e**-a - mp.mpf(1)/4
B[1] = mp.findroot(g1, A[1]*(1 + mp.mpf(10)**-15))
# 2
A[2] = mp.findroot(lambda a: (1 - mp.sech(mp.sqrt(a)))/a - mp.mpf(1)/4, mp.mpf("1.5"))
def g2(a):
    U = mp.odefun(lambda x, w: [w[1], a*w[0] - 1], 0, [mp.mpf(0), mp.mpf(0)])(1)
    V = mp.odefun(lambda x, w: [w[1], a*w[0]], 0, [mp.mpf(0), mp.mpf(1)])(1)
    c = -U[1]/V[1]
    return U[0] + c*V[0] - mp.mpf(1)/4
B[2] = mp.findroot(g2, A[2]*(1 + mp.mpf(10)**-15))
# 3
th = 4*mp.acosh(mp.e**(mp.mpf(1)/4))
A[3] = th**2*mp.e**(-mp.mpf(1)/2)/2
def g3(l):
    s = mp.findroot(lambda s: mp.odefun(lambda x, w: [w[1], -l*mp.e**w[0]], 0, [mp.mpf(0), s])(1)[0], mp.mpf("1.0"))
    return mp.odefun(lambda x, w: [w[1], -l*mp.e**w[0]], 0, [mp.mpf(0), s])(mp.mpf(1)/2)[0] - mp.mpf(1)/2
if "b3" in sys.argv:          # slow (nested odefun shooting); method B used: hd14_ref.wl
    B[3] = mp.findroot(g3, A[3]*(1 + mp.mpf(10)**-15))
# 4
Ein = lambda a: mp.ei(a) - mp.euler - mp.log(a)
A[4] = mp.findroot(lambda a: Ein(a)/(4*a) - mp.mpf(1)/2, mp.mpf("2.5"))
EinS = lambda a: mp.nsum(lambda k: a**k/(k*mp.factorial(k)), [1, mp.inf])
B[4] = mp.findroot(lambda a: EinS(a)/(4*a) - mp.mpf(1)/2, mp.mpf("2.5"))
# 5
def tsum(z):
    """Sum_odd tanh(n pi z/2)/n^5 = (1 - 2^-5) zeta(5) - Sum_odd 2/(n^5 (e^{n pi z} + 1))."""
    s = 0
    for n in range(1, 10**6, 2):
        t = 2/(mp.mpf(n)**5*(mp.e**(n*pi*z) + 1))
        s += t
        if t < mp.mpf(10)**-80: break
    return (1 - mp.mpf(2)**-5)*mp.zeta(5) - s
I5a = lambda s: s/12 - 16/pi**5*tsum(s)
I5b = lambda s: s**3/12 - 16*s**4/pi**5*tsum(1/s)
A[5] = mp.findroot(lambda s: I5a(s) - mp.mpf(1)/20, mp.mpf("1.2"))
B[5] = mp.findroot(lambda s: I5b(s) - mp.mpf(1)/20, mp.mpf("1.2"))
# 6
I6a = lambda b: mp.nsum(lambda j: 8/((2*j + 1)**4*pi**4)*(1 - mp.sech((2*j + 1)*pi*b)), [0, mp.inf])
def I6b(b):
    def term(k):
        mu = (k + mp.mpf(1)/2)*pi/b
        c = 2/(b*mu)
        return (-1)**int(k)*c/mu**2*(1 - 2*mp.tanh(mu/2)/mu)
    return mp.nsum(term, [0, mp.inf])
A[6] = mp.findroot(lambda b: I6a(b) - mp.mpf(1)/20, mp.mpf("0.3"))
B[6] = mp.findroot(lambda b: I6b(b) - mp.mpf(1)/20, mp.mpf("0.3"))
# 7
A[7] = (mp.mpf(320)/(1000*mp.sqrt(3)))**(mp.mpf(1)/4)
# 8
f8 = lambda r: 2*pi*r*((1 - r*r)/4 + r*r/2*mp.log(r)) - mp.mpf(1)/5
A[8] = mp.findroot(f8, mp.mpf("0.7"))
def g8(r):   # u'' + u'/x = -1 on [r, 1], u'(r) = 0: u = c + U(x) with U(r)=0; c = -U(1)
    U = mp.odefun(lambda x, w: [w[1], -1 - w[1]/x], r, [mp.mpf(0), mp.mpf(0)])(1)[0]
    return 2*pi*r*(-U) - mp.mpf(1)/5
B[8] = mp.findroot(g8, A[8]*(1 + mp.mpf(10)**-15))
# 9
def I9a(a):
    def term(j):
        m = 2*j + 1; k = mp.sqrt(m*m*pi*pi + a)
        return 8/(m*m*pi*pi*k*k)*(1 - 2*mp.tanh(k/2)/k)
    return mp.nsum(term, [0, mp.inf])
def I9b(a):
    inner = lambda m: mp.nsum(lambda j: 1/((2*j + 1)**2*(pi*pi*(m*m + (2*j + 1)**2) + a)), [0, mp.inf], method='richardson')
    return mp.nsum(lambda i: 64/(pi**4*(2*i + 1)**2)*inner(2*i + 1), [0, mp.inf], method='richardson')
A[9] = mp.findroot(lambda a: I9a(a) - mp.mpf(1)/50, mp.mpf("10"))
if "fast" not in sys.argv:
    mp.mp.dps = 62
    B[9] = mp.findroot(lambda a: I9b(a) - mp.mpf(1)/50, A[9])
    mp.mp.dps = 70
# 10
def I10a(a):
    def term(j):
        n = 2*j + 1; k = n*pi; D = a*a - k*k
        # g'' - k^2 g = -e^{ax}, g(0) = 0, g'(1) = 0, evaluated at x = 1 without the
        # sinh^2 - cosh^2 cancellation: g(1) = -e^a/D + a e^a tanh(k)/(D k) + 1/(D cosh k)
        g1 = -mp.e**a/D + a*mp.e**a*mp.tanh(k)/(D*k) + 1/(D*mp.cosh(k))
        return 8/(n*n*pi*pi)*g1
    return mp.nsum(term, [0, mp.inf])
def I10b(a):
    def term(j):
        mu = (j + mp.mpf(1)/2)*pi
        c = 2*(a*(-1)**int(j)*mp.e**a + mu)/(a*a + mu*mu)
        return (-1)**int(j)*c/mu**2*(1 - 2*mp.tanh(mu/2)/mu)
    return mp.nsum(term, [0, mp.inf])
A[10] = mp.findroot(lambda a: I10a(a)*mp.e**-a - mp.mpf(1)/10, mp.mpf("-0.8"))
B[10] = mp.findroot(lambda a: I10b(a)*mp.e**-a - mp.mpf(1)/10, mp.mpf("-0.8"))
for k in range(1, 11):
    s = "hd14.%d %s" % (k, mp.nstr(A[k], 66))
    if k in B: s += "  agree: %d" % agree(A[k], B[k])
    print(s, flush=True)
