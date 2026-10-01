# A non-smooth objective: a sum of absolute values whose minimum is a whole segment.
g = Abs[x - 1] + Abs[y + 2] + Abs[x - y];
FindMinimum[g, {{x, 3}, {y, 3}}]
FindMinimum[g, {{x, 3}, {y, 3}}, Method -> "NelderMead"]
