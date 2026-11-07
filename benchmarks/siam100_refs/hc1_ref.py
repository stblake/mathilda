"""hc1 (holdout A, family C1) reference, method A: mpmath at dps 80.
Each x-domain integral over (0, 1] (or (0, 2]) is mapped to a half-line
oscillatory integral (or summed from a closed form) and evaluated with
mpmath.quadosc (zeros given explicitly where the phase is t e^t: t = W(n pi)).
Method B is hc1_ref.wl (WL NIntegrate ExtrapolatingOscillatory on a
DIFFERENT substitution, or WL closed forms)."""
import mpmath as m
m.mp.dps = 80
pi = m.pi
W = lambda u: m.lambertw(u).real
res = {}
# 1  int_0^1 sin(1/x) dx = sin 1 - Ci 1
res[1] = m.sin(1) - m.ci(1)
# 2  int_0^1 sin(1/x)/x dx = pi/2 - Si(1)
res[2] = pi/2 - m.si(1)
# 3  int_0^1 cos(x^-2)/x dx = -Ci(1)/2
res[3] = -m.ci(1)/2
# 4  int_0^1 x^-1/2 cos(1/x) dx = int_1^oo t^-3/2 cos t dt
res[4] = m.quadosc(lambda t: t**(-1.5)*m.cos(t), [1, m.inf], omega=1)
# 5  int_0^1 x^-1 sin(x^-1 log x) dx = - int_0^oo sin(t e^t) dt  (x = e^-t)
res[5] = -m.quadosc(lambda t: m.sin(t*m.exp(t)), [0, m.inf], zeros=lambda n: W(n*pi))
# 6  int_0^1 cos(log x / x) dx = int_0^oo e^-t cos(t e^t) dt
res[6] = m.quadosc(lambda t: m.exp(-t)*m.cos(t*m.exp(t)), [0, m.inf],
                   zeros=lambda n: W((n - m.mpf(1)/2)*pi))
# 7  int_0^1 cos(e^(1/x)) dx = int_e^oo cos u / (u log^2 u) du  (u = e^(1/x))
res[7] = m.quadosc(lambda u: m.cos(u)/(u*m.log(u)**2), [m.e, m.inf], omega=1)
# 8  int_0^1 sin(1/(x(1-x))) dx = 2 int_4^oo sin t / (t^2 sqrt(1 - 4/t)) dt
#    and t = 4 + s^2:  = 4 int_0^oo sin(4 + s^2) / (4 + s^2)^(3/2) ds,
#    zeros s = sqrt(n pi - 4), n >= 2
g8 = lambda s: 4*m.sin(4 + s*s)/(4 + s*s)**1.5
z8 = m.sqrt(2*pi - 4)
res[8] = m.quad(g8, [0, z8]) + m.quadosc(g8, [z8, m.inf], zeros=lambda n: m.sqrt((n + 2)*pi - 4))
# 9  int_0^2 sin(1/(2-x))/(2-x) dx = int_1/2^oo sin t / t dt = pi/2 - Si(1/2)
res[9] = pi/2 - m.si(m.mpf(1)/2)
# 10 int_0^1 x^-3/2 cos(x^-1 log x) dx = int_0^oo e^(t/2) cos(t e^t) dt
res[10] = m.quadosc(lambda t: m.exp(t/2)*m.cos(t*m.exp(t)), [0, m.inf],
                    zeros=lambda n: W((n - m.mpf(1)/2)*pi))
for k in sorted(res):
    print(k, m.nstr(res[k], 70))
