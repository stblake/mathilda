# Rule chains compose automatically under fixed-point evaluation.
# The two definitions are setup lines (';' suppresses their output).
f[x_] := g[x];
g[x_] := x^2;
f[3]
