# Figure: the BFGS iterates on Rosenbrock's banana valley, collected with StepMonitor
# and drawn over the level curves of log(1 + f).
rosen = (1 - x)^2 + 100 (y - x^2)^2;
path = Prepend[Last[Reap[FindMinimum[rosen, {{x, -1.2}, {y, 1}}, StepMonitor :> Sow[{x, y}]]]][[1]], {-1.2, 1}];
Length[path]
cp = ContourPlot[Log[1 + rosen], {x, -1.6, 1.6}, {y, -0.6, 1.6}, Contours -> 25];
fig = Graphics[{First[cp], Red, Line[path], PointSize[0.008], Point[path]}, Frame -> True]
