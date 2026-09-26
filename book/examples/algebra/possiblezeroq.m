# 4.2 PossibleZeroQ: the heuristic zero test
PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]
PossibleZeroQ[Sqrt[3 + 2 Sqrt[2]] - 1 - Sqrt[2]]
PossibleZeroQ[Log[6] - Log[2] - Log[3]]
PossibleZeroQ[E^(I Pi) + 1]
PossibleZeroQ[Sin[x]^2 + Cos[x]^2]
PossibleZeroQ[x]
PossibleZeroQ[Sqrt[x^2] - x]
PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x > 0]
