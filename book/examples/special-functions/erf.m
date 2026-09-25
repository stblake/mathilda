# 4.7.3 The error functions
# Exact special values and odd symmetry
Erf[0]
Erf[Infinity]
Erf[-x]
# Arbitrary-precision real numerics track the requested precision
N[Erf[1], 20]
# The complement vanishes at infinity (cancellation-free, unlike 1 - Erf)
Erfc[Infinity]
# The imaginary error function has its own hand-rolled kernel
N[Erfi[1], 20]
# The derivative is the Gaussian, from which the Taylor series follows
D[Erf[z], z]
# The inverse error function: exact endpoints, Newton-polished numerics
InverseErf[1]
N[InverseErf[1/2], 20]
