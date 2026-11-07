"""Holdout A, family C9 (experiment 204): references for hc9.1 .. hc9.10.

Each case is the argmax (or argmin) a* of an objective F(a) built from a
parametric integral; a* is a root of F'(a).  Two independent methods each:

  hc9.1  F = (2 + cos 10a) J2(a),  J_L(a) = Int_0^L x^a sin(a/(L-x)) dx, a in [0,5]
  hc9.2  F = (2 + sin 10a) J1(a)
  hc9.10 F = (2 + sin 10a) J3(a)
         A: quadrature -- t = 1/(L-x), tanh-sinh on [1/L, 1/L+1] + quadosc on
            [1/L+1, inf) (zeros n pi/a), F' by differentiating under the integral.
         B: the Meijer-G closed form (book Sec. 9.6, generalised by x = L y)
            Int_0^1 y^a sin(b/(1-y)) dy
              = sqrt(pi) Gamma(a+1) 2^-(a+1) G^{3,0}_{2,4}(b^2/4 | (a+2)/2,(a+3)/2; 1/2,1/2,1,0),
            J_L(a) = L^(a+1) * that with b = a/L; F' by mpmath.diff.
  hc9.9  F = (2 + sin 10a) K2(a), K2(a) = Int_0^2 x^a cos(a/(2-x)) dx
         A: as above (t-substitution, quadosc), derivative under the integral.
         B: in x, the integral split at the zeros x_n = 2 - a/((n+1/2) pi) of
            cos(a/(2-x)) -- Gauss-Legendre on [0, x_0] and on each [x_n, x_{n+1}],
            the alternating tail summed by mpmath.nsum (Levin) -- derivative
            under the integral; root by findroot.
  hc9.3  argmin_{a in [0,2]} Int_0^inf x^a e^-x dx  (= Gamma(a+1))
         A: quadrature of Int x^a log(x) e^-x (= Gamma'(a+1)), findroot.
         B: digamma(a+1) = 0 (mpmath.digamma, findroot): a = 0.46163...
  hc9.4  argmax_{a in [0,4]} (1/pi) Int_0^pi cos(t - a sin t) dt  (= J1(a))
         A: quadrature of the a-derivative, findroot.  B: besseljzero(1,1,derivative=1).
  hc9.5  argmax_{a in [0,4]} (2 + sin 3a) Int_0^1 x^(-1/2) cos(a log x) dx
         (the integral = (1/2)/(1/4 + a^2))
         A: quadrature (x = e^-u: Int_0^inf e^-u/2 cos(a u) du, quadosc) of F'.  B: closed form.
  hc9.6  argmax_{b in [0,3]} Int_0^b sin(x^2) dx
         A: findroot of the integrand sin(b^2) (the derivative) near 1.77.  B: sqrt(pi).
  hc9.7  argmax_{a in [0,3]} Int_0^inf x^a e^-x sin x dx
         (= Gamma(a+1) 2^-(a+1)/2 sin((a+1) pi/4))
         A: quadrature of F' (x^a log x e^-x sin x).  B: closed form, mpmath.diff.
  hc9.8  argmin_{a in [0,10]} Int_0^inf e^-x cos(a sqrt x) dx  (= 1 - a D(a/2),
         D = Dawson's integral)
         A: quadrature of F' (-Int e^-x sqrt x sin(a sqrt x)).  B: closed form, mpmath.diff.
Usage: python3.11 hc9_ref.py hc9.1=0.64 hc9.2=0.2 hc9.10=1.42 hc9.9=0.12 hc9.5=0.15 hc9.7=1.5 hc9.8=3
"""
from mpmath import (mp, mpf, sin, cos, log, exp, sqrt, pi, gamma, meijerg, quad, quadosc, findroot,
                    diff, nstr, inf, digamma, besseljzero, atan, erfi, nsum, quadgl)

mp.dps = 75
TOL = mpf(10)**-70


def osc(g, a, L):
    """tanh-sinh on [1/L, 1/L + 1] (the x^a endpoint singularity sits at t = 1/L;
    quadosc started there loses ~7 digits), quadosc on the rest."""
    t0 = 1/mpf(L)
    return quad(g, [t0, t0 + 1]) + quadosc(g, [t0 + 1, inf], zeros=lambda n: n*pi/a)


def JL(a, L):      # Int_0^L x^a sin(a/(L-x)) dx with t = 1/(L-x): x = L - 1/t, dx = dt/t^2
    return osc(lambda t: (L - 1/t)**a*sin(a*t)/t**2, a, L)


def JLp(a, L):
    return osc(lambda t: (L - 1/t)**a*(log(L - 1/t)*sin(a*t) + t*cos(a*t))/t**2, a, L)


def KL(a, L):
    return osc(lambda t: (L - 1/t)**a*cos(a*t)/t**2, a, L)


def KLp(a, L):
    return osc(lambda t: (L - 1/t)**a*(log(L - 1/t)*cos(a*t) - t*sin(a*t))/t**2, a, L)


def Jm(a, L):
    b = a/L
    return L**(a+1)*sqrt(pi)*gamma(a+1)/2**(a+1)*meijerg([[], [(a+2)/2, (a+3)/2]],
                                                         [[mpf(1)/2, mpf(1)/2, 1], [0]], b**2/4)


def K2_x(a, deriv):
    """Int_0^2 x^a cos(a/(2-x)) dx (or its a-derivative) split at the zeros of cos(a/(2-x))."""
    if deriv:
        g = lambda x: x**a*(log(x)*cos(a/(2-x)) - sin(a/(2-x))/(2-x))
    else:
        g = lambda x: x**a*cos(a/(2-x))
    z = lambda n: 2 - a/((n + mpf(1)/2)*pi)
    n0 = 0
    while z(n0) <= 0:
        n0 += 1
    head = quad(g, [0, z(n0)])
    tail = nsum(lambda n: quad(g, [z(n0 + int(n)), z(n0 + int(n) + 1)]), [0, inf])
    return head + tail


res = []
# hc9.1, hc9.2, hc9.10: (weight, w', L, start)
for lab, wf, wpf, L, s in [("hc9.1", lambda a: 2 + cos(10*a), lambda a: -10*sin(10*a), 2, None),
                           ("hc9.2", lambda a: 2 + sin(10*a), lambda a: 10*cos(10*a), 1, None),
                           ("hc9.10", lambda a: 2 + sin(10*a), lambda a: 10*cos(10*a), 3, None)]:
    res.append((lab, wf, wpf, L))

import sys
starts = dict(a.split('=') for a in sys.argv[1:])   # e.g. hc9.1=0.6 (or hc9.2=0.2:0.21, a secant pair)
import os
ONLY = os.environ.get("ONLY", "")
want = lambda lab: not ONLY or lab in ONLY.split(",")


def st(lab):
    v = starts[lab]
    return tuple(mpf(x) for x in v.split(':')) if ':' in v else mpf(v)

out = []
for lab, wf, wpf, L in res:
    if not want(lab):
        continue
    s = st(lab)
    rA = findroot(lambda a: wpf(a)*JL(a, L) + wf(a)*JLp(a, L), s, tol=TOL)
    rB = findroot(lambda a: diff(lambda b: wf(b)*Jm(b, L), a), s, tol=TOL)
    out.append((lab, rA, rB))
    print(lab, nstr(rA, 66), nstr(rB, 66), "reldiff", nstr(abs(rA-rB)/rA, 3), flush=True)

if want("hc9.3"):
    rA = findroot(lambda a: quad(lambda x: x**a*log(x)*exp(-x), [0, 1, inf]), mpf('0.46'))
    rB = findroot(lambda a: digamma(a + 1), mpf('0.46'))
    print("hc9.3", nstr(rA, 66), nstr(rB, 66), "reldiff", nstr(abs(rA-rB)/rA, 3), flush=True)

if want("hc9.4"):
    rA = findroot(lambda a: quad(lambda t: sin(t - a*sin(t))*sin(t), [0, pi]), mpf('1.84'))
    rB = besseljzero(1, 1, derivative=1)
    print("hc9.4", nstr(rA, 66), nstr(rB, 66), "reldiff", nstr(abs(rA-rB)/rA, 3), flush=True)

if want("hc9.5"):
    I5 = lambda a: quadosc(lambda u: exp(-u/2)*cos(a*u), [0, inf], omega=a)
    I5p = lambda a: quadosc(lambda u: -u*exp(-u/2)*sin(a*u), [0, inf], omega=a)
    rA = findroot(lambda a: 3*cos(3*a)*I5(a) + (2 + sin(3*a))*I5p(a), st("hc9.5"))
    rB = findroot(lambda a: diff(lambda b: (2 + sin(3*b))*(mpf(1)/2)/(mpf(1)/4 + b*b), a), st("hc9.5"))
    print("hc9.5", nstr(rA, 66), nstr(rB, 66), "reldiff", nstr(abs(rA-rB)/rA, 3), flush=True)

if want("hc9.6"):
    rA = findroot(lambda b: sin(b*b), mpf('1.77'))
    print("hc9.6", nstr(rA, 66), nstr(sqrt(pi), 66), "reldiff", nstr(abs(rA-sqrt(pi))/rA, 3), flush=True)

if want("hc9.7"):
    rA = findroot(lambda a: quad(lambda x: x**a*log(x)*exp(-x)*sin(x), [0, 1, 5, 20, 60, inf]), st("hc9.7"))
    rB = findroot(lambda a: diff(lambda b: gamma(b+1)*2**(-(b+1)/2)*sin((b+1)*pi/4), a), st("hc9.7"))
    print("hc9.7", nstr(rA, 66), nstr(rB, 66), "reldiff", nstr(abs(rA-rB)/rA, 3), flush=True)

if want("hc9.8"):
    dawson = lambda x: sqrt(pi)/2*exp(-x*x)*erfi(x)
    rA = findroot(lambda a: quad(lambda x: -exp(-x)*sqrt(x)*sin(a*sqrt(x)), [0, 1, 5, 20, 60, inf]), st("hc9.8"))
    rB = findroot(lambda a: diff(lambda b: 1 - b*dawson(b/2), a), st("hc9.8"))
    print("hc9.8", nstr(rA, 66), nstr(rB, 66), "reldiff", nstr(abs(rA-rB)/rA, 3), flush=True)

if want("hc9.9"):
    s = st("hc9.9")
    rA = findroot(lambda a: 10*cos(10*a)*KL(a, 2) + (2 + sin(10*a))*KLp(a, 2), s, tol=TOL)
    # B is expensive (an nsum of quadratures per evaluation): start its secant at
    # A's root, so it only has to CONFIRM it (a wrong A would not be a root of B)
    print("  hc9.9 A done", nstr(rA, 30), flush=True)
    rB = findroot(lambda a: 10*cos(10*a)*K2_x(a, False) + (2 + sin(10*a))*K2_x(a, True),
                  (rA*(1 - mpf(10)**-25), rA*(1 + mpf(10)**-25)), tol=TOL)
    print("hc9.9", nstr(rA, 66), nstr(rB, 66), "reldiff", nstr(abs(rA-rB)/rA, 3), flush=True)
