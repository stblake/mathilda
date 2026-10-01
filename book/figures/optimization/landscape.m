# Figure: a one-dimensional landscape with several local minima. FindMinimum from
# x = 2, 7 and 14 descends into three different valleys (red points).
mins = FindMinimum[x Cos[x], {x, #}] & /@ {2, 7, 14}
pts = {x /. Last[#], First[#]} & /@ mins;
fig = Graphics[{First[Plot[x Cos[x], {x, 0, 16}]], Red, PointSize[0.012], Point[pts]}, Frame -> True]
