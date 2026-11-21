# NIntegrate: hard endpoints, Fourier half-lines, vertex singularities
NIntegrate[Cos[Exp[1/x]], {x, 0, 1}]
NIntegrate[Sin[x^3], {x, 0, Infinity}]
N[Gamma[1/3]/6]
NIntegrate[Exp[-x] BesselJ[0, x]^2, {x, 0, Infinity}]
NIntegrate[1/Sqrt[x^2 + y^2], {x, 0, 1}, {y, 0, 1}, PrecisionGoal -> 12]
N[2 ArcSinh[1], 13]
