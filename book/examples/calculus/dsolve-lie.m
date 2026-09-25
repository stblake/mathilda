# 4.3.7 DSolve -- the Lie point-symmetry engine
DSolve[y'[x] == 2 x y[x]/(x^2 + 2 y[x]^4 + 2), y[x], x]
DSolve[y'[x] == Sqrt[a x + b y[x] + c], y[x], x]
