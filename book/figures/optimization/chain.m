# Featured: a hanging chain of 10 links of length 1/4 between (0, 0) and (2, 0). The
# chain settles where its potential energy, the sum of the joint heights, is least.
n = 10; L = 1/4; xs = Table[Symbol["p" <> ToString[i]], {i, n - 1}]; ys = Table[Symbol["q" <> ToString[i]], {i, n - 1}];
px = Join[{0}, xs, {2}]; py = Join[{0}, ys, {0}];
links = And @@ Table[(px[[i + 1]] - px[[i]])^2 + (py[[i + 1]] - py[[i]])^2 == L^2, {i, n}];
sol = FindMinimum[Evaluate[{Total[ys], links}], Evaluate[Transpose[{Join[xs, ys], Join[Range[n - 1]/5, -0.3 + 0 Range[n - 1]]}]]];
First[sol]
fig = Graphics[{Blue, Thick, Line[Transpose[{px, py}] /. Last[sol]], Red, PointSize[0.012], Point[Transpose[{px, py}] /. Last[sol]]}, Frame -> True, AspectRatio -> Automatic]
