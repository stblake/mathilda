# Section 4.5.11: FindMinimum / FindMaximum -- local optimization from a start point.
FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1}, {y, 1}}]
FindMinimum[{x Cos[x], 1 <= x <= 15}, {x, 7}]
FindMinimum[(x - Pi)^2, {x, 0}, WorkingPrecision -> 50]
FindMaximum[Sin[x] Exp[-x/10], {x, 1}]
