# The staple operations on a square matrix, all exact. Working over the rationals,
# Inverse returns an exact rational matrix and A . Inverse[A] is the identity with
# no rounding error. A is defined once and reused across the session.
A = {{2, 1}, {1, 2}}
Det[A]
Tr[A]
Inverse[A]
A . Inverse[A]
Transpose[{{1, 2, 3}, {4, 5, 6}}]
IdentityMatrix[3]
DiagonalMatrix[{1, 2, 3}]
