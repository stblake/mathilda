# Jacobi polynomials P_n^(a,b)(x), weight (1 - x)^a (1 + x)^b
JacobiP[1, a, b, x]
JacobiP[3, 1/2, 1/2, x]
JacobiP[n, a, b, 1]
JacobiP[4, 0, 0, x] == LegendreP[4, x]
Integrate[(1 - x) (1 + x)^2 JacobiP[2, 1, 2, x] JacobiP[3, 1, 2, x], {x, -1, 1}]
D[JacobiP[n, a, b, x], x]
# Numbers come from the three-term recurrence, never the expanded polynomial
JacobiP[300, 3, 1/2, 0.9998]
N[JacobiP[300, 3, 1/2, 4999/5000], 20]
JacobiP[2.5, 1, 2, 0.3]
