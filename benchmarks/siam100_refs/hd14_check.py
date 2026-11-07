"""Holdout A, family d14: independent residual checks for hd14.3 and hd14.9
(the full second root-find is too slow at high precision).
 hd14.3: integrate the IVP u'' = -l e^u, u(0) = 0, u'(0) = th tanh(th/4) by
         mpmath.odefun (Taylor) at l = the closed-form root, dps 70; check
         u(1) = 0 and u(1/2) = 1/2.
 hd14.9: evaluate the double eigen-series (nested nsum, Richardson) at the
         closed-form root a, dps 62; check Int u - 1/50 = 0."""
import mpmath as mp
mp.mp.dps = 70
th = 4*mp.acosh(mp.e**(mp.mpf(1)/4)); l = th**2*mp.e**(-mp.mpf(1)/2)/2
s = mp.odefun(lambda x, w: [w[1], -l*mp.e**w[0]], 0, [mp.mpf(0), th*mp.tanh(th/4)])
print("hd14.3 l =", mp.nstr(l, 66), " u(1) =", mp.nstr(s(1)[0], 5), " u(1/2) - 1/2 =", mp.nstr(s(mp.mpf(1)/2)[0] - mp.mpf(1)/2, 5))
a = mp.mpf("16.0329282511843766774335528684491191374692570904225087704523220117")
mp.mp.dps = 62
pi = mp.pi
inner = lambda m: mp.nsum(lambda j: 1/((2*j + 1)**2*(pi*pi*(m*m + (2*j + 1)**2) + a)), [0, mp.inf], method='richardson')
v = mp.nsum(lambda i: 64/(pi**4*(2*i + 1)**2)*inner(2*i + 1), [0, mp.inf], method='richardson')
print("hd14.9 Int u - 1/50 at a:", mp.nstr(v - mp.mpf(1)/50, 5), "at dps 62")
