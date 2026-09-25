# 4.2 Working in a number field
# ToNumberField adjoins a generator and represents an element in the power basis.
a = ToNumberField[Sqrt[2]]
a^2
(1 + a)(1 - a)
1/(1 + a)
ToNumberField[Sqrt[2] + Sqrt[3]]
AlgebraicNumberNorm[a]
AlgebraicNumberTrace[a]
NumberFieldIntegralBasis[Sqrt[5]]
AlgebraicIntegerQ[(1 + Sqrt[5])/2]
AlgebraicIntegerQ[1/2]
