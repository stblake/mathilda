"""Holdout B, family C9 (experiment 234): references for bc9.1 .. bc9.10.
Each problem: argmax (or argmin) over a parameter interval of a parametric
integral.  Two independent methods per problem:
  A  a closed form of the integral (or of its stationarity condition) in
     special functions, root by mpmath.findroot;
  B  direct quadrature of d/da (integrand) (differentiation under the
     integral; tanh-sinh / quadosc, split at the oscillatory end), findroot.
Plus a coarse scan of F over the whole interval (printed) confirming that the
stationary point found is the GLOBAL optimum.
Run: python3.11 bc9_ref.py [k ...]   (dps 70)
"""
import sys
from mpmath import (mp, mpf, sin, cos, pi, sqrt, exp, log, gamma, ci, quad, quadosc, inf, findroot,
                    nstr, diff, angerj, besseljzero, besselj, fresnels, airyaizero, meijerg, besselyzero, bessely, linspace,
                    log10, mpc, re, erfi)
mp.dps = 70


def show(name, a, b):
    print(name, nstr(a, 66))
    print(' ' * len(name), nstr(b, 66), ' agree to', nstr(-log10(abs(a - b) / abs(a)), 3), 'digits', flush=True)


def scan(name, F, lo, hi, n=400):
    with mp.workdps(15):
        xs = linspace(mpf(lo), mpf(hi), n)
        vs = [F(x) for x in xs]
        i = max(range(n), key=lambda j: vs[j])
        print('  scan', name, 'best grid point', nstr(xs[i], 8), nstr(vs[i], 10), flush=True)


def R(f, x0):
    x0 = tuple(mpf(v) for v in x0) if isinstance(x0, tuple) else mpf(x0)
    x = findroot(f, x0, tol=mpf(10) ** -66, verify=False)
    print('   residual', nstr(abs(f(x)), 3), flush=True)
    return x
want = [int(a) for a in sys.argv[1:]] or list(range(1, 11))

if 1 in want:   # F(a) = Int_0^1 sin(a/x) dx = sin a - a Ci(a); F' = -Ci(a)
    A = R(ci, '0.6165')
    dF = lambda a: quadosc(lambda u: cos(a * u) / u, [1, inf], omega=a)   # F'(a) = Int_1^oo cos(a u)/u du
    B = R(dF, '0.6165')
    scan('bc9.1', lambda a: sin(a) - a * ci(a) if a > 0 else mpf(0), 0, 4)
    show('bc9.1', A, B)

if 2 in want:   # F(a) = (2 + sin 10a) Int_0^1 x^a sin(a/(1-x)) dx
    # A: Meijer-G closed form, K(al,be) = Int_0^1 y^al sin(be/(1-y)) dy
    #    = 2^-(al+1) sqrt(pi) Gamma(al+1) G^{3,0}_{2,4}(be^2/4 | (al+2)/2,(al+3)/2 ; 1/2,1/2,1,0)
    #    (the book's Sec. 9.6 form for Int_0^2 x^al sin(al/(2-x)) dx = 2^(al+1) K(al, al/2), extended)
    K = lambda al, be: 2 ** (-(al + 1)) * sqrt(pi) * gamma(al + 1) * meijerg(
        [[], [(al + 2) / 2, (al + 3) / 2]], [[mpf(1) / 2, mpf(1) / 2, 1], [0]], be ** 2 / 4)
    FA = lambda a: (2 + sin(10 * a)) * K(a, a)
    # B: t = 1/(1-x): Int_1^oo (1-1/t)^a sin(a t)/t^2 dt  (tanh-sinh on [1,2], quadosc on [2,oo))
    J = lambda a: quad(lambda t: (1 - 1 / t) ** a * sin(a * t) / t ** 2, [1, 2]) + \
        quadosc(lambda t: (1 - 1 / t) ** a * sin(a * t) / t ** 2, [2, inf], omega=a)
    print('  check K vs quad at a=0.8:', nstr(K(mpf('0.8'), mpf('0.8')) - J(mpf('0.8')), 5))
    FB = lambda a: (2 + sin(10 * a)) * J(a)
    with mp.workdps(20):
        scan('bc9.2', FA, '0.01', 5, 250)
    A = R(lambda a: diff(FA, a), ('0.2', '0.215'))
    B = R(lambda a: diff(FB, a), A)
    show('bc9.2', A, B)

if 3 in want:   # F(a) = (2 + cos 7a) Int_0^oo x^-a sin x dx = (2 + cos 7a) Gamma(1-a) cos(pi a/2)
    FA = lambda a: (2 + cos(7 * a)) * gamma(1 - a) * cos(pi * a / 2)
    scan('bc9.3', FA, '0.1', '1.5', 300)
    # B: quadrature; d/da x^-a sin x = -log x x^-a sin x
    I0 = lambda a: quad(lambda x: x ** -a * sin(x), [0, pi]) + quadosc(lambda x: x ** -a * sin(x), [pi, inf], omega=1)
    I1 = lambda a: quad(lambda x: -log(x) * x ** -a * sin(x), [0, pi]) + quadosc(lambda x: -log(x) * x ** -a * sin(x), [pi, inf], omega=1)
    A = R(lambda a: diff(FA, a), '0.93')
    B = R(lambda a: -7 * sin(7 * a) * I0(a) + (2 + cos(7 * a)) * I1(a), A)
    show('bc9.3', A, B)

if 4 in want:   # Anger function J_a(10) = (1/pi) Int_0^pi cos(a t - 10 sin t) dt
    FA = lambda a: angerj(a, 10)
    scan('bc9.4', FA, 0, 5, 500)
    dB = lambda a: -quad(lambda t: t * sin(a * t - 10 * sin(t)), linspace(0, pi, 9)) / pi
    A = R(lambda a: diff(FA, a), '2.14')
    B = R(dB, A)
    show('bc9.4', A, B)

if 5 in want:   # F(a) = Int_0^oo exp(-x^2) sin(a x) dx = Dawson(a/2) ; Dawson'(b) = 1 - 2 b D(b)
    D = lambda b: sqrt(pi) / 2 * exp(-b * b) * erfi(b)
    A = 2 * R(lambda b: 1 - 2 * b * D(b), '0.924')
    dB = lambda a: quad(lambda x: x * exp(-x * x) * cos(a * x), [0, 2, 4, 8, 16, inf])
    B = R(dB, '1.848')
    scan('bc9.5', lambda a: D(a / 2), 0, 10)
    show('bc9.5', A, B)

if 6 in want:   # F(a) = Int_0^oo cos(a x)/(1 + x^6) dx: residues; ARGMIN
    zs = [exp(mpc(0, 1) * pi * (2 * k + 1) / 6) for k in range(3)]
    FA = lambda a: re(pi * mpc(0, 1) * sum(exp(mpc(0, 1) * a * z) / (6 * z ** 5) for z in zs))
    print('  check residue form vs quad at a=2:', nstr(FA(mpf(2)) - quadosc(lambda x: cos(2 * x) / (1 + x ** 6), [0, inf], omega=2), 5))
    scan('bc9.6 (neg, argmin)', lambda a: -FA(a), 0, 6)
    dB = lambda a: -quadosc(lambda x: x * sin(a * x) / (1 + x ** 6), [0, inf], omega=a)
    A = R(lambda a: diff(FA, a), '3.6')
    B = R(dB, A)
    show('bc9.6', A, B)

if 7 in want:   # F(x) = Int_0^1 cos(x t)/sqrt(1-t^2) dt = (pi/2) J0(x): argmin on [0,10] = j_{1,1}
    A = besseljzero(1, 1)
    dB = lambda x: -quad(lambda th: sin(th) * sin(x * sin(th)), [0, pi / 2])   # t = sin th
    B = R(dB, '3.83')
    scan('bc7 (neg)', lambda x: -besselj(0, x), 0, 10)
    show('bc9.7', A, B)

if 8 in want:   # F(a) = Int_0^1 sin(a x^2) dx = sqrt(pi/(2a)) S(sqrt(2a/pi))
    FA = lambda a: sqrt(pi / (2 * a)) * fresnels(sqrt(2 * a / pi)) if a > 0 else mpf(0)
    scan('bc9.8', FA, 0, 20)
    A = R(lambda a: diff(FA, a), '2.2')
    dB = lambda a: quad(lambda x: x * x * cos(a * x * x), [0, 0.5, 1], method='gauss-legendre')
    B = R(dB, A)
    show('bc9.8', A, B)

if 9 in want:   # Ai(x) = (1/pi) Int_0^oo cos(t^3/3 + x t) dt; argmax on [-4,0] = first zero of Ai'
    A = airyaizero(1, derivative=1)
    # B: Ai'(-z) = (z/3) (J_{2/3}(zeta) - J_{-2/3}(zeta)), zeta = (2/3) z^(3/2)
    B = -R(lambda z: besselj(mpf(2) / 3, 2 * z ** 1.5 / 3) - besselj(-mpf(2) / 3, 2 * z ** 1.5 / 3), '1.0188')
    from mpmath import airyai
    scan('bc9.9', airyai, -4, 0)
    show('bc9.9', A, B)

if 10 in want:  # F(b) = Int_0^oo sin(x + b/x) dx = -pi sqrt(b) Y1(2 sqrt b) (checked below);
    # F' = -(pi/2) d/db [z Y1(z)], z = 2 sqrt b, and (z Y1)' = z Y0: argmax on [1/20, 3] is b = y_{0,1}^2/4
    A = besselyzero(0, 1) ** 2 / 4
    # B: F'(b) = Int_0^oo cos(x + b/x)/x dx = Int_1^oo cos(x + b/x)/x dx + Int_1^oo cos(1/u + b u)/u du
    dB = lambda b: quadosc(lambda x: cos(x + b / x) / x, [1, inf], omega=1) + \
        quadosc(lambda u: cos(1 / u + b * u) / u, [1, inf], omega=b)
    B = R(dB, '0.1997')
    FA = lambda b: -pi * sqrt(b) * bessely(1, 2 * sqrt(b))
    with mp.workdps(20):
        print('  check closed form vs quad at b=1.4:', nstr(FA(mpf('1.4')) - (
            quadosc(lambda x: sin(x + mpf('1.4') / x), [1, inf], omega=1) +
            quadosc(lambda u: sin(1 / u + mpf('1.4') * u) / u ** 2, [1, inf], omega=mpf('1.4'))), 5))
    scan('bc9.10', FA, '0.05', 3)
    show('bc9.10', A, B)
