# Section 4.5.2: the shared WorkingPrecision / AccuracyGoal / PrecisionGoal contract.
NIntegrate[Exp[-x^2], {x, 0, 1}]
NIntegrate[Exp[-x^2], {x, 0, 1}, WorkingPrecision -> 40]
NSum[1/n^2, {n, 1, Infinity}, WorkingPrecision -> 40]
