# GroebnerBasis engines: F4, multi-modular lifting, inexact and parametric bases.
cyclic5 = {a + b + c + d + e, a b + b c + c d + d e + e a, a b c + b c d + c d e + d e a + e a b, a b c d + b c d e + c d e a + d e a b + e a b c, a b c d e - 1};
Length[GroebnerBasis[cyclic5, {a, b, c, d, e}]]
GroebnerBasis[cyclic5, {a, b, c, d, e}, Method -> "F4"] === GroebnerBasis[cyclic5, {a, b, c, d, e}, Method -> "Buchberger"]
Length[GroebnerBasis[cyclic5, {a, b, c, d, e}, Modulus -> 32003, MonomialOrder -> DegreeReverseLexicographic]]
GroebnerBasis[{x^2 + y^2 - 1, x - 2 y}, {x, y}, CoefficientDomain -> InexactNumbers[20]]
GroebnerBasis[{x^2 + y^2 - r^2, x - t y}, {x, y}, CoefficientDomain -> RationalFunctions]
