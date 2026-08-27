# Section 4.5.5: ND -- numerical differentiation by Richardson extrapolation.
ND[Exp[x], x, 1]
ND[Sin[x^2], {x, 3}, 1, Terms -> 20, WorkingPrecision -> 40]
ND[Re[Cos[I y]], y, 1, AccuracyGoal -> 10]
