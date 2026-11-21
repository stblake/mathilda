# The modified Struve function, and the Anger and Weber functions
StruveL[1/2, z]
StruveL[-1/2, z]
StruveL[1, 2.]
N[StruveL[1/3, 7], 30]
# For large x, L_0 - I_0 decays like -2/(Pi x)
{StruveL[0, 30.] - BesselI[0, 30.], -2/(30 Pi) // N}
D[StruveL[nu, z], z]
# Anger and Weber: integer orders reduce to Bessel and Struve
AngerJ[2, z]
WeberE[1, z]
AngerJ[nu, 0]
{AngerJ[1/3, 2.], WeberE[1/3, 2.]}
N[AngerJ[1/3, 2], 25]
# The defining integral, by quadrature
N[1/Pi NIntegrate[Cos[t/3 - 2 Sin[t]], {t, 0, Pi}, WorkingPrecision -> 30], 25]
