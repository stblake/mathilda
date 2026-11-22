# 4.3.6 Complex line integration and the Cauchy principal value
Integrate[1/x, {x, 1 - I, 2 + 3 I}]
Integrate[z^2, {z, 0, 1 + I}]
Chop[N[Integrate[1/z, {z, 1, I, -1, -I, 1}]]]
Integrate[1/(x - 2), {x, 0, 3}, PrincipalValue -> True]
# Principal values with symbolic bounds, on the whole line, and on the half-line (Mellin)
Integrate[Exp[x]/x, {x, -a, a}, PrincipalValue -> True, Assumptions -> a > 0]
Integrate[Sin[x]/(x - 1), {x, -Infinity, Infinity}, PrincipalValue -> True]
Integrate[Exp[I a x]/(x - b), {x, -Infinity, Infinity}, PrincipalValue -> True, Assumptions -> a > 0 && Element[b, Reals]]
Integrate[x^(c - 1)/(1 - x), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> 0 < c < 1]
Integrate[Sqrt[x]/(x^2 - a^2), {x, 0, Infinity}, PrincipalValue -> True, Assumptions -> a > 0]
# Symbolic interior poles, ordered by the assumptions (v0.6683)
Integrate[1/(x - b), {x, -1, 1}, PrincipalValue -> True, Assumptions -> -1 < b < 1]
pv2 = Integrate[1/((x - a) (x - b)), {x, -1, 1}, PrincipalValue -> True, Assumptions -> -1 < a < b < 1];
Simplify[pv2]
N[{pv2 /. {a -> -1/3, b -> 1/2}, Integrate[1/((x + 1/3) (x - 1/2)), {x, -1, 1}, PrincipalValue -> True]}]
