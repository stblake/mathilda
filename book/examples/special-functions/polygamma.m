# 4.7.1 LogGamma and the polygamma family
# LogGamma keeps integers exact as a log of a factorial
LogGamma[10]
# ... and stays finite exactly where Gamma overflows a machine double
Gamma[171.]
LogGamma[171.]
# The digamma function at a positive integer: a rational minus Euler's constant
PolyGamma[5]
# Gauss's digamma theorem closes rational arguments into elementary form
PolyGamma[0, 1/2]
# The trigamma at 1 is Zeta[2]
PolyGamma[1, 1]
# LogGamma differentiates to digamma; each derivative raises the polygamma order
D[LogGamma[z], z]
D[PolyGamma[0, z], z]
