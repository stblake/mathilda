# DSolve: three newer families -- Riemann P, y' given implicitly, scaling symmetry
DSolve[x (1 - x) y''[x] + (1 - 3 x) y'[x] - y[x] == 0, y[x], x]
DSolve[Sin[y'[x]] + y'[x] == x, y[x], x]
DSolve[4 y'[x]^3 + 4 y'[x] == x, y[x], x]
DSolve[x y[x] y''[x] - x y'[x]^2 + y[x] y'[x] == 0, y[x], x]
