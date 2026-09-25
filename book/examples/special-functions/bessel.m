# 4.7.4 Bessel functions of the first kind and its siblings
# J_0(0) = 1; J_n(0) = 0 for n > 0
BesselJ[0, 0]
# Half-integer orders collapse to elementary trigonometric form
BesselJ[1/2, x]
BesselJ[-1/2, x]
# Arbitrary-precision numerics
N[BesselJ[0, 1], 20]
# The Frobenius series about the origin
Series[BesselJ[0, x], {x, 0, 6}]
# The recurrence surfaces in the derivative
D[BesselJ[0, x], x]
# The three companions: second kind Y, and modified I and K (machine precision)
BesselY[0, 1.]
BesselI[0, 1.]
BesselK[0, 1.]
