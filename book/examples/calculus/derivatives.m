# 4.3.1 Derivatives -- D, Derivative, Dt
# The unknown-function examples come first, while f and g are still free.
D[x^3, x]
D[Sin[x^2], x]
D[Sin[a x], {x, 3}]
D[Log[b, x], x]
D[f[g[x]], x]
Derivative[2][f][x]
D[x^2 + 5 y^3, {{x, y}}]
D[x^2 + 5 y^3, {{x, y}, 2}]
D[{x^2 + y, x y}, {{x, y}}]
Dt[y^2 + Sin[x]]
Dt[x^n, x]
# Now give f a definition, and its derivative operator computes.
f[x_] := x^5 + 6 x^3
f'[x]
f'[5]
