# 4.4.5 Rank, row reduction, and the null space
RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
MatrixRank[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
NullSpace[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
m . First[NullSpace[m]]
RowReduce[{{a, b, c}, {d, e, f}, {a + d, b + e, c + f}}]
MatrixRank[{{a, b}, {2 a, 2 b}}]
