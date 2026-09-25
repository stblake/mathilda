# 4.7.5 Legendre functions of the first and second kind
# The Legendre polynomials P_n, generated on demand
Table[LegendreP[n, x], {n, 0, 3}]
LegendreP[4, x]
# Orthogonality on [-1, 1]: the squared norm is 2/(2n+1)
Integrate[LegendreP[2, x]^2, {x, -1, 1}]
# ... and distinct degrees are orthogonal
Integrate[LegendreP[2, x] LegendreP[3, x], {x, -1, 1}]
# The second-kind solution LegendreQ carries the logarithm
LegendreQ[0, x]
LegendreQ[1, x]
LegendreQ[2, x]
# It evaluates numerically like any other
N[LegendreQ[0, 1/2], 20]
