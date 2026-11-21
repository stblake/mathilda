# The Kelvin functions: Bessel functions on the rotated ray x E^(3 Pi I/4)
KelvinBer[z]
{KelvinBer[0, 0], KelvinBei[0, 0], KelvinKei[0, 0]}
KelvinBer[2.] + I KelvinBei[2.]
BesselJ[0, 2. Exp[3 Pi I/4]]
N[KelvinKer[1/3, -2 + 3 I], 20]
D[KelvinBer[z], z]
