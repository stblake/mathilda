"""C9 reference: argmax_{a in [0,5]} (2+sin 10a) J(a), J(a) = Int_0^2 x^a sin(a/(2-x)) dx.
Substitute t = 1/(2-x): J = Int_{1/2}^inf (2-1/t)^a sin(a t) t^-2 dt, split as
tanh-sinh quad on [1/2, 1] (the x^a endpoint singularity sits at t = 1/2; quadosc
started there loses ~7 digits) + quadosc on [1, inf) with the zeros n pi/a.
J' by differentiating under the integral sign; root of F'(a) by findroot (secant)."""
from mpmath import mp, mpf, sin, cos, log, quad, quadosc, findroot, inf, nstr, pi
mp.dps = 75
def osc(g, a): return quad(g, [mpf(1)/2, 1]) + quadosc(g, [1, inf], zeros=lambda n: n*pi/a)
def J(a):  return osc(lambda t: (2-1/t)**a*sin(a*t)/t**2, a)
def Jp(a): return osc(lambda t: (2-1/t)**a*(log(2-1/t)*sin(a*t) + t*cos(a*t))/t**2, a)
Fp = lambda a: 10*cos(10*a)*J(a) + (2+sin(10*a))*Jp(a)
r = findroot(Fp, mpf('0.7859336743'), tol=mpf(10)**-130)
print(nstr(r, 68))
