# Hankel functions and the spherical Bessel family
HankelH1[2.5, 3.]
BesselJ[2.5, 3.] + I BesselY[2.5, 3.]
HankelH1[0, 0]
{SphericalBesselJ[0, 0], SphericalBesselJ[2, 0]}
SphericalBesselJ[2, z]
N[SphericalBesselJ[2, 3/2], 30]
# The elementary form j_2(z) = (3/z^3 - 1/z) Sin[z] - 3 Cos[z]/z^2
N[(3/(3/2)^3 - 1/(3/2)) Sin[3/2] - 3/(3/2)^2 Cos[3/2], 30]
{SphericalHankelH1[1, 2.], SphericalHankelH2[1, 2.]}
D[SphericalBesselJ[n, z], z]
# FunctionExpand writes the half-integer family in sines and cosines
FunctionExpand[SphericalBesselJ[2, z]]
FunctionExpand[{SphericalBesselY[1, z], SphericalHankelH1[0, z]}]
FunctionExpand[HankelH1[1/2, z]]
