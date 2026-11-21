# 4.4.9 Arbitrary-precision dense linear algebra
MatrixRank[N[HilbertMatrix[40]]]
MatrixRank[N[HilbertMatrix[40], 100]]
LinearSolve[N[HilbertMatrix[12], 40], N[ConstantArray[1, 12], 40]][[1]]
Eigenvalues[N[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}, 30]]
Eigenvalues[{N[{{2, 1}, {1, 3}}, 30], N[{{1, 0}, {0, 2}}, 30]}]
SingularValueList[N[HilbertMatrix[4], 30]]
