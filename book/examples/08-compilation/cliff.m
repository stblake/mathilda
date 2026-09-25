# Section 8.4: the compilable subset is a cliff. One unlowerable head anywhere
# in the body sends the WHOLE body back to the interpreter.
CompileDiagnostics[{{x, _Real}}, Sin[x] + Gamma[x]]
CompileDiagnostics[{{x, _Real}}, Sin[x] + BarnesG[x]]
CompileDiagnostics[{{v, _Real, 1}}, Total[Riffle[v, 0.]]]
CompileDiagnostics[{{n, _Integer}}, MoebiusMu[n] + EulerPhi[n]]
