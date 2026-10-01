# Figure: minimize the distance to (2, 1) over the unit disk. The optimum is where a level
# curve of the objective first touches the circle.
sol = FindMinimum[{(x - 2)^2 + (y - 1)^2, x^2 + y^2 <= 1}, {x, y}, Method -> "SLSQP"]
cp = ContourPlot[(x - 2)^2 + (y - 1)^2, {x, -1.5, 2.5}, {y, -1.5, 2}, Contours -> 15];
fig = Graphics[{First[cp], Black, Thick, Circle[{0, 0}, 1], Red, PointSize[0.014], Point[{x, y} /. Last[sol]], Point[{2, 1}]}, Frame -> True, AspectRatio -> Automatic]
