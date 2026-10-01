# Figure: Himmelblau's function has four minima. FindMinimum reaches each one from a
# start in the matching quadrant.
h = (x^2 + y - 11)^2 + (x + y^2 - 7)^2;
mins = Table[{x, y} /. Last[FindMinimum[h, {{x, s[[1]]}, {y, s[[2]]}}]], {s, {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}}}]
cp = ContourPlot[Log[1 + h], {x, -5, 5}, {y, -5, 5}, Contours -> 20];
fig = Graphics[{First[cp], Red, PointSize[0.012], Point[mins]}, Frame -> True]
