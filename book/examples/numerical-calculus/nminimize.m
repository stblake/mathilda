# Section 4.5.11: NMinimize / NMaximize -- global optimization, constraints, integers.
NMinimize[20 + x^2 - 10 Cos[2 Pi x] + y^2 - 10 Cos[2 Pi y], {x, y}]
NMinimize[{x + y, x^2 + y^2 <= 1}, {x, y}]
NMinimize[{x + y, x + 2 y >= 3 && x >= 0 && y >= 0, Element[{x, y}, Integers]}, {x, y}]
NMinimize[{20 + x^2 - 10 Cos[2 Pi x] + y^2 - 10 Cos[2 Pi y], -5.12 <= x <= 5.12 && -5.12 <= y <= 5.12}, {x, y}, Method -> "DIRECT"]
