# Holdout-A D6 family (experiment 211) references, method A:
# mpmath tanh-sinh with the interval split at every zero of the oscillating
# factor (the phase = k*pi or (k+1/2)*pi), dps 80; cross-checked with
# Gauss-Legendre on the same split.  Method B is hd6_ref.wl.
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi
def split(a, b, phase_inv, ks):
    pts = [mp.mpf(a)]
    for k in ks:
        x = phase_inv(k)
        if a < x < b: pts.append(x)
    pts.append(mp.mpf(b)); return sorted(pts)
P = {}
# 1  cos(exp(x^2)) on [0,2]
P[1] = (lambda x: mp.cos(mp.exp(x**2)), split(0, 2, lambda k: mp.sqrt(mp.log((k+mp.mpf(1)/2)*pi)), range(0, 20)))
# 2  exp(x + cos(exp(exp(x)))) on [-1,1]
P[2] = (lambda x: mp.exp(x + mp.cos(mp.exp(mp.exp(x)))), split(-1, 1, lambda k: mp.log(mp.log((k+mp.mpf(1)/2)*pi)), range(0, 8)))
# 3  sin(exp(exp(2x))) on [0,1]
P[3] = (lambda x: mp.sin(mp.exp(mp.exp(2*x))), split(0, 1, lambda k: mp.log(mp.log(k*pi))/2, range(1, 520)))
# 4  sin(x^3) on [0,3]
P[4] = (lambda x: mp.sin(x**3), split(0, 3, lambda k: mp.cbrt(k*pi), range(1, 10)))
# 5  exp(sin(exp(x))) on [0,pi]
P[5] = (lambda x: mp.exp(mp.sin(mp.exp(x))), split(0, pi, lambda k: mp.log(k*pi), range(1, 9)))
# 6  sin(exp(4x))/(2+x) on [-1,1]
P[6] = (lambda x: mp.sin(mp.exp(4*x))/(2+x), split(-1, 1, lambda k: mp.log(k*pi)/4, range(1, 20)))
# 7  x exp(x^2) sin(exp(x^2)) on [0,2]
P[7] = (lambda x: x*mp.exp(x**2)*mp.sin(mp.exp(x**2)), split(0, 2, lambda k: mp.sqrt(mp.log(k*pi)), range(1, 20)))
# 8  exp(x + sin(exp(exp(x+1/2)))) on [-1,1]
P[8] = (lambda x: mp.exp(x + mp.sin(mp.exp(mp.exp(x + mp.mpf(1)/2)))), split(-1, 1, lambda k: mp.log(mp.log(k*pi)) - mp.mpf(1)/2, range(1, 30)))
# 9  cos(x^2 + x) on [0,8]
P[9] = (lambda x: mp.cos(x**2 + x), split(0, 8, lambda k: (-1 + mp.sqrt(1 + 4*(k+mp.mpf(1)/2)*pi))/2, range(0, 25)))
# 10 cos(1/(x+1/100)) on [0,1]
a = mp.mpf(1)/100
P[10] = (lambda x: mp.cos(1/(x+a)), split(0, 1, lambda k: 1/((k+mp.mpf(1)/2)*pi) - a, range(0, 40)))
closed = {
  7: (mp.cos(1) - mp.cos(mp.exp(4)))/2,
  10: (lambda F: F(1/a) - F(1/(1+a)))(lambda u: -mp.cos(u)/u - mp.si(u)),
  4: mp.quad(lambda x: mp.sin(x**3), [0, 3]) * 0,   # placeholder, closed form in hd6_ref.wl
}
for k in range(1, 11):
    f, pts = P[k]
    A = mp.quad(f, pts)
    B = mp.quad(f, pts, method='gauss-legendre')
    line = "hd6.%d  %s  TS-GL %s" % (k, mp.nstr(A, 70), mp.nstr(A - B, 3))
    if k in (7, 10):
        line += "  closed-diff %s" % mp.nstr(A - closed[k], 3)
    print(line, flush=True)
