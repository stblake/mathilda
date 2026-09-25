# 4.4.8 Eigenvalues, eigenvectors, and the Jordan form
CharacteristicPolynomial[{{1, 2}, {3, 4}}, x]
Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Eigenvalues[{{a, b}, {c, d}}]
Eigenvalues[{{0, -1}, {1, 0}}]
Eigenvectors[{{2, 1, 0}, {0, 2, 0}, {0, 0, 3}}]
m = {{27, 48, 81}, {-6, 0, 0}, {1, 0, 3}};
{sj, jj} = JordanDecomposition[m];
jj
m == sj . jj . Inverse[sj]
