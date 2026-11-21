# Weierstrass elliptic functions with invariants {g2, g3}
WeierstrassP[u, {0, 0}]
# A double root of 4t^3 - g2 t - g3 degenerates to a hyperbolic function
WeierstrassP[u, {12, -8}]
WeierstrassP[0.3, {1., 2.}]
N[WeierstrassP[3/10, {1, 2}], 30]
# The differential equation wp'^2 = 4 wp^3 - g2 wp - g3, at 30 digits
p = WeierstrassP[3/10`30, {1, 2}]; pp = WeierstrassPPrime[3/10`30, {1, 2}];
pp^2 - (4 p^3 - p - 2)
# The period lattice, and double periodicity
{w1, w3} = WeierstrassHalfPeriods[{2`30, 3`30}]
Abs[{WeierstrassP[3/10 + 2 w1, {2, 3}], WeierstrassP[3/10 + 2 w3, {2, 3}]} - WeierstrassP[3/10, {2, 3}]]
WeierstrassInvariants[{w1, w3}]
# Derivatives close on the family
D[WeierstrassPPrime[u, {g2, g3}], u]
D[WeierstrassZeta[u, {g2, g3}], u]
D[WeierstrassSigma[u, {g2, g3}], u]
