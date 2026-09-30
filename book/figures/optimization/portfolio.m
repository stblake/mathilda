# Featured: the efficient frontier of three assets. For each target return, minimize the
# portfolio variance w.C.w over long-only weights that sum to 1.
mu = {0.06, 0.10, 0.14}; cov = {{0.010, 0.002, 0.001}, {0.002, 0.040, 0.010}, {0.001, 0.010, 0.090}};
risk[r_] := Sqrt[First[FindMinimum[{{a, b, c} . cov . {a, b, c}, a + b + c == 1 && mu . {a, b, c} == r && a >= 0 && b >= 0 && c >= 0}, {{a, 0.3}, {b, 0.3}, {c, 0.4}}, Method -> "SLSQP"]]];
frontier = Table[{risk[r], r}, {r, 0.065, 0.135, 0.005}];
fig = ListPlot[frontier, Joined -> True]
