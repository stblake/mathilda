# Special functions compile at arbitrary precision too
f = Compile[{{x, _Real}}, Erf[x] x^2 + BesselJ[0, x], WorkingPrecision -> 30];
f[1.5`30]
N[Erf[3/2] 9/4 + BesselJ[0, 3/2], 30]
CompileDiagnostics[{{x, _Real}}, MathieuC[2, 1, x], WorkingPrecision -> 30]
