# DSolve: the PDE initial-value layer (v0.6496)
DSolve[{D[u[x, t], t] == D[u[x, t], {x, 2}], u[x, 0] == E^(-x^2)}, u[x, t], {x, t}]
DSolve[{D[u[x, t], {t, 2}] == c^2 D[u[x, t], {x, 2}], u[x, 0] == E^(-x^2), Derivative[0, 1][u][x, 0] == 0}, u[x, t], {x, t}]
DSolve[{D[u[x, t], t] + u[x, t] D[u[x, t], x] == 0, u[x, 0] == x}, u[x, t], {x, t}]
DSolve[{D[u[x, t], t] + u[x, t] D[u[x, t], x] + D[u[x, t], {x, 3}] == 0, u[x, 0] == 3 Sech[x/2]^2}, u[x, t], {x, t}]
DSolve[{D[u[y, t], t] == nu D[u[y, t], {y, 2}], u[y, 0] == 0, u[0, t] == U}, u[y, t], {y, t}]
DSolve[D[u[x, y], x, y] == E^u[x, y], u[x, y], {x, y}]
