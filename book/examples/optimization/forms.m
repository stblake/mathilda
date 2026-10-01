# Local versus global: the same objective, x cos x, under the four optimizers.
FindMinimum[x Cos[x], {x, 2}]
FindMinimum[x Cos[x], {x, 7}]
NMinimize[{x Cos[x], 0 <= x <= 16}, x]
FindMaximum[x Cos[x], {x, 5}]
NMaximize[{x Cos[x], 0 <= x <= 16}, x]
{fmin, rules} = FindMinimum[x Cos[x], {x, 2}]
x /. rules
x Cos[x] /. rules
