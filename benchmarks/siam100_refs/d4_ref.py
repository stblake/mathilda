# D4 reference, method B: an expansion of R(e^-u) for large u that never
# touches the Gram/Moebius series.  With R(x) = sum_{k>=0} (log x)^k /(k k! zeta(k+1))
# (k = 0 term := 1) and Ramanujan's master theorem,
#   R(e^-u) = (1/2 pi i) Int_{(c)} Gamma(s) u^-s / (-s zeta(1-s)) ds,  0 < c < 1/2;
# shifting the contour to the right collects the poles of 1/zeta(1-s):
#   s = 2n+1 (trivial zeros):  (-1)^(n+1) 2 (2 pi)^(2n) / ((2n+1) zeta(2n+1) u^(2n+1))
#   s = 1-rho (nontrivial zeros): -2 Re[ Gamma(1-rho) u^(rho-1) / ((1-rho) zeta'(rho)) ]
# The first series converges geometrically (ratio (2 pi/u)^2); the second
# like exp(-pi gamma/2), so 70 zero pairs give ~ 100 digits.  The root is
# found in t = log x = -u by mpmath.findroot (secant) at dps 110.
import mpmath as mp
mp.mp.dps = 110
NZ = 70
zz = []
for j in range(1, NZ + 1):
    rho = mp.zetazero(j)
    zz.append((rho, mp.zeta(rho, derivative=1)))
def R(u):
    s = mp.mpf(0)
    n = 1
    while True:
        term = (-1)**(n+1) * 2 * (2*mp.pi)**(2*n) / ((2*n+1) * mp.zeta(2*n+1) * u**(2*n+1))
        s += term
        if abs(term) < mp.mpf(10)**(-mp.mp.dps - 20): break
        n += 1
    for rho, zp in zz:
        s -= 2 * mp.re(mp.gamma(1-rho) * u**(rho-1) / ((1-rho) * zp))
    return s
print("check vs WL u^3 R(e^-u) at u=34000 (0.32838933605...):", mp.nstr(R(mp.mpf(34000))*34000**3, 15))
print("last zero term size:", mp.nstr(abs(mp.gamma(1-zz[-1][0]) / zz[-1][1]), 5))
u = mp.findroot(R, (mp.mpf(34100), mp.mpf(34150)), solver='secant', tol=mp.mpf(10)**-200)
print("t* =", mp.nstr(-u, 80))
print("x* =", mp.nstr(mp.exp(-u), 75))
