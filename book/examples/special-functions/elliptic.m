# 4.7.6 Legendre elliptic integrals
# At m = 0 the integrand is 1, so K is the quarter period itself
EllipticK[0]
# A lemniscatic singular value, closed through the gamma function
EllipticK[1/2]
# The same number to thirty digits
N[EllipticK[1/2], 30]
# The integrand loses its last denominator at t = Pi/2
EllipticK[1]
# The parameter expansion at m = 0, from a dedicated Series kernel
Series[EllipticK[m], {m, 0, 2}]
# A full quarter period closes the incomplete form into the complete one
EllipticF[Pi/2, 1/2]
# The amplitude derivative of the incomplete integral IS its integrand
D[EllipticF[phi, m], phi]
# The quasi-period: F(phi + k Pi | m) = F(phi | m) + 2 k K(m), as a residual
N[EllipticF[1/5 + 3 Pi, 1/2] - (EllipticF[1/5, 1/2] + 6 EllipticK[1/2]), 20]
# The parameter derivative of the second kind, in terms of E and K themselves
D[EllipticE[m], m]
# E(phi | 1) is the integral of Abs[Cos[t]], which is Sin[phi] on |phi| <= Pi/2
EllipticE[1/2, 1]
# Beyond Pi/2 that identity is false, and the true value is 2 - Sin[2]
N[EllipticE[2, 1] - (2 - Sin[2]), 20]
# The third kind reduces to the first when its characteristic vanishes
EllipticPi[0, phi, m]
# Past the pole at Sin[t]^2 == 1/n the value is complex, not a principal value
N[EllipticPi[3/2, 1/2], 20]
# The perimeter of an ellipse with a = 1 and eccentricity e^2 = 1/4
N[4 EllipticE[1/4], 20]
# A pendulum released at 30 degrees swings 1.74% slower than the small-angle law
N[2 EllipticK[Sin[Pi/12]^2]/Pi, 20]
# K is a Gauss hypergeometric function in disguise, to 25 digits
N[Pi/2 Hypergeometric2F1[1/2, 1/2, 1, 1/2] - EllipticK[1/2], 25]
