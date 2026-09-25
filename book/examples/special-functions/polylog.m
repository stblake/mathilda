# 4.7.6 Polylogarithm, Lerch transcendent, and the Lambert W function
# The polylogarithm generalizes the logarithm (order 1) ...
PolyLog[1, x]
# ... and the dilogarithm has famous closed forms
PolyLog[2, 1]
PolyLog[2, -1]
PolyLog[2, 1/2]
# At argument 1 it is the zeta function
PolyLog[3, 1]
# The Lerch transcendent contains all of these as special cases
LerchPhi[z, 1, 1]
# The Lambert W function (ProductLog) inverts w e^w = z
ProductLog[E]
ProductLog[-1/E]
N[ProductLog[1], 20]
# Its derivative is written implicitly in terms of itself
D[ProductLog[x], x]
