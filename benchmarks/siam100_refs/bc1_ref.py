"""Holdout B, family C1 (experiment 226) -- reference values, method A (mpmath).

For each of the ten endpoint-/infinite-range oscillatory integrals:
  * a closed form where one exists (special functions in mpmath), and
  * mpmath.quadosc on the integral in a substituted form (t = 1/x for the
    x -> 0 problems), summing the pieces between the zeros by nsum.
Both at dps 80; the two must agree to >= 60 digits.  bc1.3 and bc1.10 have
no closed form: their second method is PARI/gp (siam100_refs/bc1_ref.gp,
sumalt over Gauss pieces between zeros), compared by hand.
"""
import mpmath as mp

mp.mp.dps = 80
pi = mp.pi
R = {}

def both(k, closed, osc):
    a, b = closed, osc
    d = -mp.log10(abs(a - b) / abs(a)) if a != b else 80
    print("bc1.%d" % k, mp.nstr(a, 70), " agree %.1f" % d)
    R[k] = a

# 1. Int_0^1 sin(1/x) dx = Int_1^oo sin t / t^2 dt = sin 1 - Ci(1)
both(1, mp.sin(1) - mp.ci(1),
     mp.quadosc(lambda t: mp.sin(t) / t**2, [1, mp.inf], omega=1))
# 2. Int_0^1 sin(1/x)/x dx = Int_1^oo sin t / t dt = pi/2 - Si(1)
both(2, pi / 2 - mp.si(1),
     mp.quadosc(lambda t: mp.sin(t) / t, [1, mp.inf], omega=1))
# 3. Int_0^1 cos(log(x)/x) dx = Int_0^oo cos(t e^t) e^-t dt  (x = e^-t)
#    no closed form; zeros t e^t = (n - 1/2) pi.  Second method: PARI.
v3 = mp.quadosc(lambda t: mp.cos(t * mp.exp(t)) * mp.exp(-t), [0, mp.inf],
                zeros=lambda n: mp.lambertw((n - mp.mpf(1) / 2) * pi).real)
v3b = mp.quadosc(lambda t: mp.cos(t * mp.exp(t)) * mp.exp(-t), [0, mp.inf],
                 zeros=lambda n: mp.lambertw(n * pi).real)  # other zero set (sin zeros)
both(3, v3, v3b)
# 4. Int_0^oo sin(x^2) dx = sqrt(pi/8)
both(4, mp.sqrt(pi / 8),
     mp.quadosc(lambda x: mp.sin(x**2), [0, mp.inf], zeros=lambda n: mp.sqrt(n * pi)))
# 5. Int_0^oo cos(x)/sqrt(x) dx = sqrt(pi/2)
both(5, mp.sqrt(pi / 2),
     mp.quadosc(lambda u: 2 * mp.cos(u**2), [0, mp.inf],
                zeros=lambda n: mp.sqrt((n - mp.mpf(1) / 2) * pi)))
# 6. Int_0^1 sin(1/x)/sqrt(x) dx = Int_1^oo sin t t^(-3/2) dt
#    = 2 sin 1 + 2 Int_1^oo cos t / sqrt t dt
#    = 2 sin 1 + 2 (sqrt(pi/2) - 2 sqrt(pi/2) C(sqrt(2/pi)))   (C: Fresnel, pi t^2/2)
c6 = 2 * mp.sin(1) + 2 * (mp.sqrt(pi / 2) - 2 * mp.sqrt(pi / 2) * mp.fresnelc(mp.sqrt(2 / pi)))
both(6, c6, mp.quadosc(lambda t: mp.sin(t) * t**(-mp.mpf(3) / 2), [1, mp.inf], omega=1))
# 7. Int_0^oo sin(x)/(x+1) dx = Ci(1) sin 1 + (pi/2 - Si(1)) cos 1
both(7, mp.ci(1) * mp.sin(1) + (pi / 2 - mp.si(1)) * mp.cos(1),
     mp.quadosc(lambda x: mp.sin(x) / (x + 1), [0, mp.inf], omega=1))
# 8. Int_0^oo J1(x)/sqrt(x) dx = 2^(-1/2) Gamma(3/4)/Gamma(5/4)
#    (Weber: Int x^mu J_nu = 2^mu Gamma((nu+mu+1)/2)/Gamma((nu-mu+1)/2))
both(8, mp.gamma(mp.mpf(3) / 4) / mp.gamma(mp.mpf(5) / 4) / mp.sqrt(2),
     mp.quad(lambda u: 2 * mp.besselj(1, u**2), [0, 1])        # x = u^2 on [0, 1]
     + mp.quadosc(lambda x: mp.besselj(1, x) / mp.sqrt(x), [1, mp.inf], omega=1))
# 9. Int_0^oo cos(x^3) dx = Gamma(4/3) cos(pi/6)
both(9, mp.gamma(mp.mpf(4) / 3) * mp.cos(pi / 6),
     mp.quadosc(lambda x: mp.cos(x**3), [0, mp.inf],
                zeros=lambda n: mp.cbrt((n - mp.mpf(1) / 2) * pi)))
# 10. Int_2^oo sin(x)/log(x) dx: no closed form.  Two mpmath zero sets here,
#     PARI sumalt (bc1_ref.gp) as the independent method.
v10 = mp.quadosc(lambda x: mp.sin(x) / mp.log(x), [2, mp.inf], omega=1)
v10b = mp.quadosc(lambda x: mp.sin(x) / mp.log(x), [2, mp.inf],
                  zeros=lambda n: n * pi)
both(10, v10, v10b)
