# Higher-level mathematics ships as rules and programs in Mathilda's own language.
# HeavisideTheta's values are DownValues loaded from src/internal/distributions.m:
Length[DownValues[HeavisideTheta]]
{HeavisideTheta[3], HeavisideTheta[-2], HeavisideTheta[0]}
Context[HeavisideTheta]
# Derivatives, by contrast, moved into C (src/calculus/deriv.c):
DownValues[D]
D[f[g[x]], x]
# Integration methods such as Trager's algorithm are .m packages, loaded on first use:
Length[Names["Trager`*"]]
Integrate[x/Sqrt[1 + x^2], x, Method -> "Trager"]
Length[Names["Trager`*"]]
