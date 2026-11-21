# 4.4.5 Elimination over radical and parametric entries
m = {{Sqrt[c], c, a}, {c, c^(3/2), a Sqrt[c]}, {1, b, 0}};
MatrixRank[m]
NullSpace[m]
m . First[NullSpace[m]] // Together
LinearSolve[{{Sqrt[c], 1}, {c, Sqrt[c] + 1}}, {1, 0}]
