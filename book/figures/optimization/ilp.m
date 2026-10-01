# Figure: an integer program. The shaded triangle is the feasible region of the linear
# relaxation, whose optimum (blue) is the far vertex (4, 4.5); the integer optimum (red)
# is the lattice point (1, 2), nowhere near it.
lattice = Flatten[Table[{i, j}, {i, 0, 5}, {j, 0, 5}], 1];
fig = Graphics[{Opacity[0.3], Blue, Polygon[{{0, 0.5}, {4, 4.5}, {0, 1.3}}], Opacity[1], Gray, PointSize[0.012], Point[lattice], Blue, PointSize[0.02], Point[{4, 4.5}], Red, Point[{1, 2}]}, Frame -> True, AspectRatio -> Automatic]
