# Hardy's Z function and the Riemann-Siegel theta function
{RiemannSiegelTheta[0], RiemannSiegelZ[0]}
RiemannSiegelTheta[10.]
N[RiemannSiegelTheta[10], 25]
N[Im[LogGamma[1/4 + 5 I]] - 5 Log[Pi], 25]
# Z is real on the real line, and |Z(t)| = |Zeta(1/2 + I t)|
RiemannSiegelZ[{1., 20.5}]
Abs[Zeta[1/2 + 20.5 I]]
# Riemann's R function against the prime-counting function
RiemannR[2]
N[RiemannR[10^6], 20]
{PrimePi[10^6], N[LogIntegral[10^6]]}
