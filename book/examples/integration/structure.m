# Symbolic Integration -- the structure theorem folds dependent kernels into the tower
Integrate[Log[x] Log[2 x], x]
f = D[ArcTan[E^x]^3, x]
Integrate[f, x]
Risch`ElementaryIntegralQ[E^Tanh[x], x]
Risch`ElementaryIntegralQ[E^x (Erf[x] + Erfc[x]), x]
