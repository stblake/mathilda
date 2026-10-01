# The minimum of the Gamma function, at machine and at 40-digit precision.
FindMinimum[Gamma[x], {x, 1.5}]
FindMinimum[Gamma[x], {x, 1.5}, WorkingPrecision -> 40]
FindRoot[PolyGamma[x], {x, 1.5}, WorkingPrecision -> 40]
