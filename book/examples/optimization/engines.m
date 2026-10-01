# The global engines side by side on Rastrigin, then over five seeds each.
rast = 20 + x^2 - 10 Cos[2 Pi x] + y^2 - 10 Cos[2 Pi y];
NMinimize[rast, {x, y}, Method -> "SimulatedAnnealing"]
NMinimize[rast, {x, y}, Method -> "DualAnnealing"]
NMinimize[rast, {x, y}, Method -> "BasinHopping"]
NMinimize[rast, {x, y}, Method -> "SHGO"]
NMinimize[rast, {x, y}, Method -> "RandomSearch"]
NMinimize[rast, {x, y}, Method -> {"RandomSearch", "SearchPoints" -> 400}]
Table[First[NMinimize[rast, {x, y}, Method -> {"RandomSearch", "RandomSeed" -> s}]], {s, 1, 5}]
Table[First[NMinimize[rast, {x, y}, Method -> {"SimulatedAnnealing", "RandomSeed" -> s}]], {s, 1, 5}]
