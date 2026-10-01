# Featured: the largest open box from a 30 x 20 sheet by cutting squares of side x
# from the corners and folding up the sides.
vol = (30 - 2 x) (20 - 2 x) x;
sol = FindMaximum[{vol, 0 <= x <= 10}, {x, 2}]
fig = Graphics[{First[Plot[vol, {x, 0, 10}]], Red, PointSize[0.014], Point[{x /. Last[sol], First[sol]}]}, Frame -> True]
