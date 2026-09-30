# Featured: the cheapest can. Minimize the metal (surface area) of a closed cylinder
# holding 355 ml. With h eliminated the area is 2 Pi r^2 + 710/r.
sol = NMinimize[{2 Pi r^2 + 2 Pi r h, Pi r^2 h == 355 && r >= 1 && h >= 1}, {r, h}]
{r0, a0} = {r /. Last[sol], First[sol]};
fig = Graphics[{First[Plot[2 Pi r^2 + 710/r, {r, 1.5, 8}]], Red, PointSize[0.014], Point[{r0, a0}]}, Frame -> True]
