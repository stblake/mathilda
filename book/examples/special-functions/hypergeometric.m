# 4.7.6 The hypergeometric family: the common ancestor of the elementary functions
# The Gauss function 2F1 reduces to a logarithm for these parameters
Hypergeometric2F1[1, 1, 2, x]
# A negative integer numerator parameter terminates the series into a polynomial
Hypergeometric2F1[-3, 1, 1, x]
# Kummer's confluent 1F1 gives (E^x - 1)/x here ...
Hypergeometric1F1[1, 2, x]
# ... and the exponential itself when the parameters agree
Hypergeometric1F1[a, a, x]
# The confluent 0F1 is a hyperbolic cosine of a square root
Hypergeometric0F1[1/2, x]
# The generalized pFq carries any number of parameters; here it matches the 2F1 above
HypergeometricPFQ[{1, 1}, {2}, x]
# A numeric value (this one is 2 Log[2])
N[Hypergeometric2F1[1, 1, 2, 1/2], 20]
