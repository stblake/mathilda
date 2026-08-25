# Solving A x = b, and the structure of the solution set: LinearSolve solves it,
# RowReduce is the Gaussian elimination behind it (reduced row echelon form), and
# MatrixRank / NullSpace report the rank and a basis for what the matrix annihilates.
LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]
RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]
MatrixRank[{{1, 2}, {2, 4}}]
NullSpace[{{1, 2}, {2, 4}}]
