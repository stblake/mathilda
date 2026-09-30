# Featured: Fermat's principle. Light runs from (0, 1) in air (speed 1) to (2, -1) in glass
# (speed 2/3) and crosses the surface y = 0 where the travel time is least.
time = Sqrt[x^2 + 1] + Sqrt[(2 - x)^2 + 1]/(2/3);
sol = FindMinimum[time, {x, 1}]
x0 = x /. Last[sol];
{Sin[ArcTan[x0]], Sin[ArcTan[2 - x0]], Sin[ArcTan[x0]]/Sin[ArcTan[2 - x0]]}
fig = Graphics[{Gray, Line[{{-0.5, 0}, {2.5, 0}}], Thick, Red, Line[{{0, 1}, {x0, 0}, {2, -1}}], Blue, Line[{{0, 1}, {2, -1}}]}, Frame -> True, AspectRatio -> Automatic]
