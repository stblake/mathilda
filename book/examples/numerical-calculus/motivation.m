# Section 4.5.1: where symbolic methods stall, numerical methods still deliver a number.
Integrate[Exp[-x^2], {x, 0, 1}]
NIntegrate[Exp[-x^2], {x, 0, 1}]
Integrate[Sin[x]/x, {x, 0, 1}]
NIntegrate[Sin[x]/x, {x, 0, 1}]
Solve[x^5 - x - 1 == 0, x]
NRoots[x^5 - x - 1 == 0, x]
FindRoot[Cos[x] == x, {x, 1}]
