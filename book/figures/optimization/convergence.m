# Figure: convergence of four local methods on Rosenbrock from (-1.2, 1). Each curve is
# log10 of the objective at successive iterates.
rosen = (1 - x)^2 + 100 (y - x^2)^2;
trace[m_] := Log[10, rosen /. Thread[{x, y} -> #]] & /@ Last[Reap[FindMinimum[rosen, {{x, -1.2}, {y, 1}}, Method -> m, StepMonitor :> Sow[{x, y}]]]][[1]];
data = trace /@ {"QuasiNewton", "Newton", "LBFGSB", "TrustExact"};
Length /@ data
fig = ListPlot[data, Joined -> True]
