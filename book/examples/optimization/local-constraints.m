# Constrained local minimization: boxes, a disk, a corner, an equality, and HS71.
FindMinimum[{(x - 3)^2 + (y - 3)^2, 0 <= x <= 1 && 0 <= y <= 2}, {x, y}]
FindMinimum[{(x - 2)^2 + (y - 1)^2, x^2 + y^2 <= 1}, {x, y}]
FindMinimum[{x + y, 3 x + 2 y >= 7 && x >= 0 && y >= 0}, {x, y}]
FindMinimum[{x + y, 3 x + 2 y >= 7 && x >= 0 && y >= 0}, {x, y}, Method -> "SLSQP"]
FindMinimum[{x y, x^2 + y^2 == 1}, {{x, 1}, {y, 0.5}}]
FindMaximum[{x y, x + 2 y <= 4 && x >= 0 && y >= 0}, {x, y}]
FindMinimum[{x1 x4 (x1 + x2 + x3) + x3, x1 x2 x3 x4 >= 25 && x1^2 + x2^2 + x3^2 + x4^2 == 40 && 1 <= x1 <= 5 && 1 <= x2 <= 5 && 1 <= x3 <= 5 && 1 <= x4 <= 5}, {{x1, 1}, {x2, 5}, {x3, 5}, {x4, 1}}, Method -> "SLSQP"]
FindMinimum[{(x - 2)^2 + (y - 1)^2, x^2 + y^2 <= 1}, {x, y}, Method -> "COBYQA"]
