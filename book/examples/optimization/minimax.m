# Linear constraints by SQP, minimax objectives, and constraints at WorkingPrecision
FindMinimum[{(x - 3)^2 + (y - 2)^2, x + y <= 1 && x - y >= -2}, {{x, 0}, {y, 0}}]
FindMinimum[Max[Abs[x - 1], Abs[y + 2], Abs[x - y]], {{x, 0}, {y, 0}}]
cheb = Max[Table[Abs[Exp[t] - (a + b t)], {t, 0, 1, 1/10}]];
FindMinimum[cheb, {{a, 1}, {b, 1}}]
N[E - 1]
NMinimize[cheb, {a, b}]
FindMinimum[{(x - 3)^2 + (y - 2)^2, x + y <= 1}, {{x, 0}, {y, 0}}, WorkingPrecision -> 30]
FindMinimum[{x^2 + 2 y^2, x^2 + y^2 >= 1}, {{x, 1}, {y, 1}}, WorkingPrecision -> 30]
FindMinimum[Sum[(x[i] - i)^2 + Cos[x[i]], {i, 1, 3}], Table[{x[i], i}, {i, 1, 3}], WorkingPrecision -> 40]
