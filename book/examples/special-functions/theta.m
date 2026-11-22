# Jacobi theta functions of nome q
EllipticTheta[3, u, 0]
EllipticTheta[1, 0, q]
EllipticTheta[1, 0.7, 0.3]
N[EllipticTheta[1, 1/2, 1/10], 30]
D[EllipticTheta[3, u, q], u]
# Jacobi's quartic identity and theta_1' = theta_2 theta_3 theta_4, at 30 digits
q = 1/10`30;
EllipticTheta[3, q]^4 - EllipticTheta[2, q]^4 - EllipticTheta[4, q]^4
EllipticThetaPrime[1, 0, q] - EllipticTheta[2, q] EllipticTheta[3, q] EllipticTheta[4, q]
