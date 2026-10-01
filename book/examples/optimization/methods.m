# The local-method menu on Rosenbrock's function, then Himmelblau's four minima.
rosen = (1 - x)^2 + 100 (y - x^2)^2;
FindMinimum[rosen, {{x, -1.2}, {y, 1}}, Method -> "Newton"]
FindMinimum[rosen, {{x, -1.2}, {y, 1}}, Method -> "LBFGSB"]
FindMinimum[rosen, {{x, -1.2}, {y, 1}}, Method -> "TrustExact"]
FindMinimum[rosen, {{x, -1.2}, {y, 1}}, Method -> "NelderMead"]
h = (x^2 + y - 11)^2 + (x + y^2 - 7)^2;
FindMinimum[h, {{x, 1}, {y, 1}}]
FindMinimum[h, {{x, -1}, {y, 1}}]
FindMinimum[h, {{x, -1}, {y, -1}}]
FindMinimum[h, {{x, 1}, {y, -1}}]
FindMaximum[h, {{x, 0}, {y, 0}}]
