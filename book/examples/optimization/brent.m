# One variable: Brent's method, the bracketing forms, and evaluation counts.
FindMinimum[x^4 - 3 x^2 + x, {x, 1}]
FindMinimum[x^4 - 3 x^2 + x, {x, -1}]
FindMinimum[x Cos[x], {x, 7, 1, 15}]
FindMinimum[Abs[x - 1/3], {x, 0}]
n = 0; FindMinimum[Cos[x] + x/5, {x, 3}, EvaluationMonitor :> n++]; n
