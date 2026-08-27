# Section 4.5.9: NRoots -- every root of a polynomial, real and complex.
NRoots[1 + 2 x + 3 x^2 + 4 x^3 == 0, x]
NRoots[x^5 - x - 1 == 0, x]
NRoots[x^2 - (3 + 4 I) == 0, x]
NRoots[(x - 1)^3 == 0, x]
NRoots[Product[x - k, {k, 1, 15}] == 0, x]
NRoots[x^2 - 2 == 0, x, PrecisionGoal -> 30]
