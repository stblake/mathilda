# Gegenbauer (ultraspherical) polynomials C_n^lambda(x)
GegenbauerC[3, a, x]
{GegenbauerC[n, 1, x], GegenbauerC[n, 1/2, x]}
GegenbauerC[n, x]
D[GegenbauerC[n, a, x], x]
GegenbauerC[2.5, 1.5, 0.3]
N[GegenbauerC[5/2, 3/2, 3/10], 30]
# The hypergeometric definition, independently
N[Gamma[5/2 + 3]/(Gamma[3] Gamma[5/2 + 1]) Hypergeometric2F1[-5/2, 5/2 + 3, 2, (1 - 3/10)/2], 30]
