# 4.3.7 DSolve -- the first-order family
DSolve[y'[x] == y[x], y[x], x]
DSolve[y'[x] == x y[x], y[x], x]
DSolve[y'[x] + y[x]/x == x^2, y[x], x]
DSolve[y'[x] == y[x] + x y[x]^2, y[x], x]
DSolve[2 x y[x] + (x^2 + 1) y'[x] == 0, y[x], x]
DSolve[y[x] == x y'[x] + y'[x]^2, y[x], x]
DSolve[y'[x] == (x + y[x])^2, y[x], x]
DSolve[{y'[x] == 3 y[x], y[0] == 5}, y[x], x]
