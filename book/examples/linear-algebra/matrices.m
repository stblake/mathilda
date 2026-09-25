# 4.4.1 Vectors and matrices as expressions
# A matrix is a List of Lists -- an ordinary expression, so every generic tool
# (Part, Map, Transpose) works on it with no special matrix type.
m = {{2, 1, 3}, {0, 5, 4}, {6, 0, 7}}
Dimensions[m]
MatrixQ[m]
m[[All, 2]]
Map[Total, m]
Transpose[m]
ConjugateTranspose[{{1, 2 + I}, {3 - I, 4 I}}]
IdentityMatrix[3]
DiagonalMatrix[{a, b, c}]
Diagonal[m]
