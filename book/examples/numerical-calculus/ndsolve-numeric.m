# Section 4.5.8: NDSolve returns InterpolatingFunctions you can evaluate anywhere.
NDSolve[{y'[x] == -1000 (y[x] - Cos[x]) - Sin[x], y[0] == 1}, y, {x, 0, 3}, Method -> "BackwardEuler"]
sol = NDSolve[{y''[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 6}]
(y[x] /. sol[[1]]) /. x -> 2.0
Cos[2.0]
NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 6}]
