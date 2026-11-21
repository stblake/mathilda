# Mathieu functions: y'' + (a - 2 q Cos[2 z]) y = 0
# At q = 0 the equation is y'' + a y = 0
MathieuC[a, 0, z]
MathieuS[a, 0, z]
# Machine and arbitrary-precision values
MathieuC[2, 1, 3.2]
N[MathieuS[2, 1, 37/10], 30]
MathieuC[2, 1, 3.2`30 + 0.5`30 I]
# Characteristic values a_r(q), b_r(q)
{MathieuCharacteristicA[2, 1.], MathieuCharacteristicB[2, 1.]}
N[MathieuCharacteristicA[0, 25], 30]
# At a_2(q) the even solution is pi-periodic and DLMF-normalised
a2 = N[MathieuCharacteristicA[2, 1], 30];
MathieuC[a2, 1, 2/5 + Pi] - MathieuC[a2, 1, 2/5]
NIntegrate[MathieuC[MathieuCharacteristicA[2, 1.], 1., z]^2, {z, 0, 2 Pi}]
# At a_2(q) no periodic odd solution exists: the Floquet-normalised odd function is 0
MathieuS[MathieuCharacteristicA[2, 1.], 1., 0.5]
# The derivative closes on the pair
D[MathieuCPrime[a, q, x], x]
# The Wronskian C S' - C' S is constant
Table[MathieuC[2., 1., x] MathieuSPrime[2., 1., x] - MathieuCPrime[2., 1., x] MathieuS[2., 1., x], {x, {0.3, 1.7, 5.2}}]
Series[MathieuC[3, 2, x], {x, 0, 6}]
