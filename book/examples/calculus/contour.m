# 4.3.6 Complex line integration and the Cauchy principal value
Integrate[1/x, {x, 1 - I, 2 + 3 I}]
Integrate[z^2, {z, 0, 1 + I}]
Chop[N[Integrate[1/z, {z, 1, I, -1, -I, 1}]]]
Integrate[1/(x - 2), {x, 0, 3}, PrincipalValue -> True]
