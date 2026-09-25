# 4.4.7 Matrix decompositions: LU, QR, SVD, Schur
LUDecomposition[{{1, 1, 1}, {2, 4, 8}, {3, 9, 27}}]
{q, r} = QRDecomposition[{{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}}];
q
r
Transpose[q] . r
SingularValueDecomposition[{{1, 2}, {1, 2}}]
schur = {{2., 1.}, {1., 3.}};
{qs, ts} = SchurDecomposition[schur];
ts
Chop[schur - qs . ts . ConjugateTranspose[qs]]
