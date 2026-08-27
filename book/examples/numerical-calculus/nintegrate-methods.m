# Section 4.5.3: choosing an NIntegrate Method for the shape of the integrand.
NIntegrate[Exp[-x^2], {x, 0, Infinity}, Method -> "DoubleExponential"]
NIntegrate[Cos[100 x] Exp[-x], {x, 0, Infinity}, Method -> "LevinRule"]
NIntegrate[1/x, {x, -1, 2}, Method -> "PrincipalValue", Exclusions -> {0}]
NIntegrate[Exp[-x^2 - y^2], {x, -3, 3}, {y, -3, 3}, Method -> "MonteCarlo"]
