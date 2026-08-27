# Section 4.5.10: FindRoot -- Newton, secant, Brent, and systems.
FindRoot[Cos[x] == x, {x, 1}]
FindRoot[Cos[x] == x, {x, 1}, WorkingPrecision -> 50]
FindRoot[Cos[x] - x, {x, 0.5, 0, 1}]
FindRoot[Sin[x] + Exp[x], {x, 0, -1}]
FindRoot[{x^2 + y^2 == 4, x y == 1}, {{x, 2}, {y, 0.5}}]
