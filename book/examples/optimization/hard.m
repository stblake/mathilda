# Harder landscapes with published optima: six-hump camel, drop-wave, Schwefel, Eggholder.
NMinimize[{(4 - 2.1 x^2 + x^4/3) x^2 + x y + (-4 + 4 y^2) y^2, -3 <= x <= 3 && -2 <= y <= 2}, {x, y}]
NMinimize[{-(1 + Cos[12 Sqrt[x^2 + y^2]])/(0.5 (x^2 + y^2) + 2), -5.12 <= x <= 5.12 && -5.12 <= y <= 5.12}, {x, y}]
NMinimize[{418.9829*2 - x Sin[Sqrt[Abs[x]]] - y Sin[Sqrt[Abs[y]]], -500 <= x <= 500 && -500 <= y <= 500}, {x, y}]
egg = -(y + 47) Sin[Sqrt[Abs[x/2 + y + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]];
NMinimize[{egg, -512 <= x <= 512 && -512 <= y <= 512}, {x, y}]
NMinimize[{egg, -512 <= x <= 512 && -512 <= y <= 512}, {x, y}, Method -> "SimulatedAnnealing"]
