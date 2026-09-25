# 4.7.3 The integral functions: where non-elementary antiderivatives get names
# Integrals with no elementary form return the special function that names them
Integrate[Sin[x]/x, x]
Integrate[Cos[x]/x, x]
Integrate[Exp[x]/x, x]
Integrate[1/Log[x], x]
# Each evaluates to full precision
N[SinIntegral[2], 20]
N[CosIntegral[2], 20]
N[ExpIntegralEi[2], 20]
N[LogIntegral[2], 20]
# Si saturates to Pi/2 at infinity (the Dirichlet integral)
SinIntegral[Infinity]
# The Fresnel integrals of optics and diffraction
N[FresnelC[1], 20]
N[FresnelS[1], 20]
# Si differentiates back to the cardinal sine
D[SinIntegral[z], z]
