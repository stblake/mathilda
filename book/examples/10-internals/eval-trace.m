# Trace shows every rewrite the evaluator makes on the way to a fixed point.
f[x_] := x^2
Trace[f[3]]
Trace[2 + 3 + 4]
