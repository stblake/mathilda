# Featured: a production plan. Two products earn 3 and 5 per batch; three plants limit
# the batches to x <= 4, 2 y <= 12 and 3 x + 2 y <= 18. The best plan is a vertex.
sol = NMaximize[{3 x + 5 y, x <= 4 && 2 y <= 12 && 3 x + 2 y <= 18 && x >= 0 && y >= 0}, {x, y}]
fig = Graphics[{Opacity[0.3], Blue, Polygon[{{0, 0}, {4, 0}, {4, 3}, {2, 6}, {0, 6}}], Opacity[1], Gray, Line[{{0, 7.2}, {6, 3.6}}], Line[{{0, 4.8}, {6, 1.2}}], Red, PointSize[0.018], Point[{x, y} /. Last[sol]]}, Frame -> True, AspectRatio -> Automatic]
