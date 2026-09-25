# 4.3.7 DSolve -- partial differential equations
DSolve[D[u[t, x], t] + c D[u[t, x], x] == 0, u, {t, x}]
DSolve[D[u[x, y], x] + 3 D[u[x, y], y] + u[x, y] == 1, u, {x, y}]
DSolve[{D[u[x, t], {t, 2}] == c^2 D[u[x, t], {x, 2}], u[x, 0] == Sin[x], Derivative[0, 1][u][x, 0] == 0}, u, {x, t}]
DSolve[{D[u[x, t], t] == k D[u[x, t], {x, 2}], u[x, 0] == f[x]}, u, {x, t}]
DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x]
