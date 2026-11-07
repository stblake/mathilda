"""C9 reference, second method (independent of c9_ref.py's quadrature):
the closed form of the book (Sec. 9.6, from wolframscript's Integrate)
  J(a) = Int_0^2 x^a sin(a/(2-x)) dx = sqrt(pi) Gamma(a+1) G^{3,0}_{2,4}(a^2/16 | (a+2)/2, (a+3)/2 ; 1/2, 1/2, 1, 0),
evaluated with mpmath.meijerg; F(a) = (2 + sin 10a) J(a); root of F'(a)
(mpmath.diff, numerical) by findroot."""
from mpmath import mp, mpf, sin, sqrt, pi, gamma, meijerg, findroot, diff, nstr
mp.dps = 90
J = lambda a: sqrt(pi)*gamma(a+1)*meijerg([[], [(a+2)/2, (a+3)/2]], [[mpf(1)/2, mpf(1)/2, 1], [0]], a**2/16)
print("J(0.785933674350371) =", nstr(J(mpf('0.785933674350371')), 25))
F = lambda a: (2 + sin(10*a))*J(a)
r = findroot(lambda a: diff(F, a), mpf('0.7859336743'), tol=mpf(10)**-150)
print(nstr(r, 70))
