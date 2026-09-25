# 4.3.7 DSolve -- variable-coefficient and special-function equations
DSolve[y''[x] - x y[x] == 0, y[x], x]
DSolve[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y[x], x]
DSolve[x^2 y''[x] + 4 x y'[x] + 7 y[x] == 0, y[x], x]
DSolve[y''[x] - (x^2 + 3/(4 x^2)) y[x] == 0, y[x], x]
DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x]
