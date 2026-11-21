# Appell functions of two variables
AppellF1[a, b1, b2, c, x, 0]
AppellF1[a, b1, b2, c, x, x]
AppellF1[a, b1, b2, a, x, y]
AppellF1[-1, 1, 2, 3, 1/2, 1/3]
AppellF1[0.3, 0.2, 0.5, 1.7, 0.4, 0.3]
N[AppellF1[3/10, 1/5, 1/2, 17/10, 2/5, 3/10], 30]
# Picard's Euler integral for F1, by quadrature
N[Gamma[17/10]/(Gamma[3/10] Gamma[7/5]) NIntegrate[t^(-7/10) (1 - t)^(2/5) (1 - 2/5 t)^(-1/5) (1 - 3/10 t)^(-1/2), {t, 0, 1}, WorkingPrecision -> 30], 25]
# Far outside the double-series region, by analytic continuation
N[AppellF1[3/10, 1/5, 1/2, 17/10, 10, -20], 20]
{AppellF2[0.5, 1., 1.5, 2., 2.5, 0.2, 0.3], AppellF3[0.5, 1., 1.5, 2., 2.5, 0.2, 0.3], AppellF4[0.5, 1., 1.5, 2., 0.04, 0.09]}
D[AppellF4[a, b, c1, c2, x, y], x]
# The incomplete elliptic integral F is an F1
{EllipticF[0.7, 0.3], Sin[0.7] AppellF1[1/2, 1/2, 1/2, 3/2, Sin[0.7]^2, 0.3 Sin[0.7]^2]}
