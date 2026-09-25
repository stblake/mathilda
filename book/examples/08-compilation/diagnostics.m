# Section 8.2: CompileDiagnostics reports the result type and instruction counts.
CompileDiagnostics[{{x, _Real}}, x^2 + 1]
CompileDiagnostics[{{x, _Real}}, Sin[x]^2 + Cos[x]^2]
CompileDiagnostics[{{v, _Real, 1}}, Total[v] + Max[v]]
