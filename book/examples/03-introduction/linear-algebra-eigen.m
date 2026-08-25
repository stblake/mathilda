# Eigenvalues and eigenvectors, computed exactly. The golden ratio falls out of a
# 2x2 matrix, whose powers are the Fibonacci numbers; and when the eigenvalues have
# no radical form they come back as Root objects, exactly as Solve does for a cubic.
# Then two ways to the same numbers: numericalise the exact Root eigenvalues with N,
# or compute the eigenvalues of the floating-point matrix directly.
Eigenvalues[{{2, 1}, {1, 2}}]
Eigenvectors[{{2, 1}, {1, 2}}]
Eigenvalues[{{1, 1}, {1, 0}}]
MatrixPower[{{1, 1}, {1, 0}}, 10]
Eigenvalues[{{0, 1, 0}, {0, 0, 1}, {1, 1, 0}}]
Eigenvalues[{{0, 1, 0}, {0, 0, 1}, {1, 1, 0}}] // N
Eigenvalues[{{0, 1, 0}, {0, 0, 1}, {1, 1, 0}} // N]
