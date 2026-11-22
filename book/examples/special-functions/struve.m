# The Struve function H_nu
StruveH[0, 0]
StruveH[1, 2.]
N[StruveH[0, 1], 30]
# For large x, H_0 - Y_0 approaches 2/(Pi x)
StruveH[0, 10.] - BesselY[0, 10.]
N[2/(10 Pi)]
# The first zero of J_0
BesselJZero[0, 1]
N[BesselJZero[0, 1], 30]
BesselJ[0, N[BesselJZero[0, 1], 30]]
# Half-integer order is elementary; the derivative follows the recurrence
StruveH[1/2, x]
StruveH[-1/2, x]
D[StruveH[nu, x], x]
# At order 0 the recurrence closes on itself; a negative integer order reduces to a positive one
D[StruveH[0, x], x]
StruveH[-2, x]
