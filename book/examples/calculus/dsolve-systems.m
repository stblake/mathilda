# 4.3.7 DSolve -- systems of ODEs
DSolve[{y'[x] == y[x] - 2 z[x], z'[x] == y[x] - z[x]}, {y[x], z[x]}, x]
DSolve[{x'[t] == x[t] + y[t], y'[t] == 4 x[t] + y[t]}, {x[t], y[t]}, t]
