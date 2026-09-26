# Quantifier elimination: Resolve / ForAll / Exists
Resolve[Exists[x, x^2 == 2], Reals]
Resolve[Exists[x, x^2 == -1], Reals]
Resolve[ForAll[x, x^2 + 1 > 0], Reals]
Resolve[Exists[x, x^2 + b x + c == 0], Reals]
Resolve[ForAll[x, Exists[y, x + y == 0]], Reals]
Resolve[ForAll[eps, eps > 0, Exists[del, del > 0, del < eps]], Reals]
