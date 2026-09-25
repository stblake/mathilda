# 4.4.3 Determinant, inverse, and matrix powers
Det[{{1, 2}, {3, 4}}]
Det[{{a, b, c}, {d, e, f}, {g, h, i}}]
Inverse[{{1, 2}, {3, 4}}]
Inverse[{{a, b}, {c, d}}]
p = {{1, 2}, {3, 4}};
p . Inverse[p]
MatrixPower[{{1, 1}, {1, 0}}, 10]
MatrixPower[{{a, b}, {c, d}}, 2]
