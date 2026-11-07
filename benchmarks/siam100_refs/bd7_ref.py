# Holdout B, family D7 (areas of planar regions via root-finding + 1-D
# integrals with algebraic endpoint/interior singularities).
# Method A (this file): mpmath findroot + tanh-sinh with the roots/kinks as
# breakpoints, dps 90.  Closed forms where they exist (bd7.1, bd7.2).
# Method B: bd7_ref.gp (PARI/gp solve on brackets + intnum with the endpoint
# singularity exponents declared, realprecision 90).
import mpmath as mp
mp.mp.dps = 90
F = mp.findroot; Q = mp.quad
cbrt = lambda t: mp.sign(t) * mp.cbrt(abs(t))
out = {}
# 1. lens x^2 < y < cos x
x0 = F(lambda x: mp.cos(x) - x**2, mp.mpf("0.8"))
out["bd7.1"] = Q(lambda x: mp.cos(x) - x**2, [-x0, x0]); cf1 = 2*(mp.sin(x0) - x0**3/3)
# 2. e^x - 2 < y < 1 - x^2
a = F(lambda x: mp.exp(x) + x**2 - 3, mp.mpf("-1.7")); b = F(lambda x: mp.exp(x) + x**2 - 3, mp.mpf("0.8"))
out["bd7.2"] = Q(lambda x: 3 - x**2 - mp.exp(x), [a, b]); cf2 = (3*b - b**3/3 - mp.exp(b)) - (3*a - a**3/3 - mp.exp(a))
# 3. sin x sin y > 1/2 in [0,pi]^2
out["bd7.3"] = Q(lambda x: mp.pi - 2*mp.asin(1/(2*mp.sin(x))), [mp.pi/6, mp.pi/2, 5*mp.pi/6])
# 4. exp(-(y - x^2)^2) > 1/2 + y^2/8
g4 = lambda y: mp.mpf(1)/2 + y**2/8; s4 = lambda y: mp.sqrt(-mp.log(g4(y)))
yA = F(lambda y: -mp.log(g4(y)) - y**2, mp.mpf("-0.75")); yB = F(lambda y: -mp.log(g4(y)) - y**2, mp.mpf("0.75"))
out["bd7.4"] = Q(lambda y: 2*mp.sqrt(y + s4(y)), [yA, 0, yB]) + Q(lambda y: 2*(mp.sqrt(y + s4(y)) - mp.sqrt(y - s4(y))), [yB, 2])
# 5. bounded oval of y^2 < x^3 - 3x + 1
r1 = F(lambda x: x**3 - 3*x + 1, mp.mpf("-1.9")); r2 = F(lambda x: x**3 - 3*x + 1, mp.mpf("0.35"))
out["bd7.5"] = Q(lambda x: 2*mp.sqrt(x**3 - 3*x + 1), [r1, 0, r2])
# 6. cos x > x^2 + y^4
out["bd7.6"] = Q(lambda x: 2*mp.root(mp.cos(x) - x**2, 4), [-x0, 0, x0])
# 7. sin(x+y) > x^2 + y^2 - 1, u = (x+y)/sqrt 2
h7 = lambda u: 1 + mp.sin(mp.sqrt(2)*u) - u**2
u1 = F(h7, mp.mpf("-0.55")); u2 = F(h7, mp.mpf("1.4"))
out["bd7.7"] = Q(lambda u: 2*mp.sqrt(h7(u)), [u1, 0, u2])
# 8. exp(-(2y + x^3)^2) > y^2/16 + exp(cos y - 1)
g8 = lambda y: y**2/16 + mp.exp(mp.cos(y) - 1); s8 = lambda y: mp.sqrt(-mp.log(g8(y)))
w8 = lambda y: cbrt(s8(y) - 2*y) - cbrt(-s8(y) - 2*y)
p = F(lambda y: g8(y) - 1, mp.mpf("3.7")); m = F(lambda y: g8(y) - 1, mp.mpf("-3.7"))
out["bd7.8"] = Q(w8, [m, 0, p])
# 9. x^2 + y^2 < 9/4 and cos x + cos y > 1
c9 = lambda x: mp.sqrt(mp.mpf(9)/4 - x**2); a9 = lambda x: mp.acos(1 - mp.cos(x))
x1 = F(lambda x: c9(x) - a9(x), mp.mpf("0.75")); x2 = F(lambda x: c9(x) - a9(x), mp.mpf("1.3"))
out["bd7.9"] = 4*(Q(c9, [0, x1]) + Q(a9, [x1, x2]) + Q(c9, [x2, mp.mpf(3)/2]))
# 10. exp(-(y + x^3)^2) > y^2/100 + exp(sin y): three islands
g10 = lambda y: y**2/100 + mp.exp(mp.sin(y)); s10 = lambda y: mp.sqrt(-mp.log(g10(y)))
w10 = lambda y: cbrt(-y + s10(y)) - cbrt(-y - s10(y))
G = lambda y: g10(y) - 1
y1 = F(G, mp.mpf("-7.94")); y2 = F(G, mp.mpf("-7.04")); y3 = F(G, mp.mpf("-3.04"))
yc = F(lambda y: -mp.log(g10(y)) - y**2, mp.mpf("-0.86"))
y4 = F(G, mp.mpf("3.25")); y5 = F(G, mp.mpf("5.85"))
out["bd7.10"] = Q(w10, [y1, y2]) + Q(w10, [y3, yc, 0]) + Q(w10, [y4, y5])
for k in sorted(out, key=lambda s: int(s.split(".")[1])):
    print(k, mp.nstr(out[k], 75))
print("bd7.1 closed", mp.nstr(cf1, 75)); print("bd7.2 closed", mp.nstr(cf2, 75))
print("# roots: x0", mp.nstr(x0, 20), "a,b", mp.nstr(a, 20), mp.nstr(b, 20), "yA,yB", mp.nstr(yA, 20), mp.nstr(yB, 20),
      "r1,r2", mp.nstr(r1, 20), mp.nstr(r2, 20), "u", mp.nstr(u1, 20), mp.nstr(u2, 20), "p,m", mp.nstr(p, 20), mp.nstr(m, 20),
      "x1,x2", mp.nstr(x1, 20), mp.nstr(x2, 20), "y1..y5,yc", [mp.nstr(t, 20) for t in (y1, y2, y3, y4, y5, yc)])
