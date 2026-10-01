# format: png
# Figure: the Eggholder function on [-512, 512]^2. Differential evolution stops in a deep
# local basin (white); multi-chain simulated annealing reaches the global minimum at the
# edge of the box (red).
egg = -(y + 47) Sin[Sqrt[Abs[x/2 + y + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]];
box = -512 <= x <= 512 && -512 <= y <= 512;
de = {x, y} /. Last[NMinimize[{egg, box}, {x, y}]];
sa = {x, y} /. Last[NMinimize[{egg, box}, {x, y}, Method -> "SimulatedAnnealing"]];
dp = DensityPlot[egg, {x, -512, 512}, {y, -512, 512}, PlotPoints -> 150];
fig = Graphics[{First[dp], White, Disk[de, 12], Red, Disk[sa, 12]}, Frame -> True]
