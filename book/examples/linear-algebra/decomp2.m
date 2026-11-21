# 4.4.7 More decompositions: Cholesky, Hessenberg, singular values, LU back-substitution
c = CholeskyDecomposition[{{4, 2}, {2, 3}}]
ConjugateTranspose[c] . c
CholeskyDecomposition[{{1, 2}, {2, 1}}]
{p, h} = HessenbergDecomposition[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 10.}}];
h
Chop[p . h . ConjugateTranspose[p] - {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 10.}}]
SingularValueList[{{1, 2}, {3, 4}}]
SingularValueList[{{1., 2.}, {2., 4.}}]
lu = LUDecomposition[{{1, 2}, {3, 4}}];
LUBackSubstitution[lu, {5, 6}]
LUBackSubstitution[lu, {{1, 0}, {0, 1}}]
