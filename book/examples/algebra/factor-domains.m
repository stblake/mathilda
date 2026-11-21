# Factor over other coefficient domains, and the sign convention.
Factor[x^4 + 1]
Factor[x^4 + 1, Modulus -> 2]
Factor[x^4 + 1, Modulus -> 3]
Factor[x^2 - 2 y^2, Extension -> Sqrt[2]]
Factor[x^2 + y^2, GaussianIntegers -> True]
Factor[E^(2 x) - 1]
Length[FactorList[x^1260 - 1]] - 1
Factor[(1 - 4 x) (1 + 3 y)]
