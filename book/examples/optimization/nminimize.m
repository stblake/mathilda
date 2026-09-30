# Global search: Rastrigin's function defeats a local method and yields to NMinimize.
rast = 20 + x^2 - 10 Cos[2 Pi x] + y^2 - 10 Cos[2 Pi y];
FindMinimum[rast, {{x, 2.3}, {y, -1.6}}]
NMinimize[rast, {x, y}]
NMinimize[{rast, -5.12 <= x <= 5.12 && -5.12 <= y <= 5.12}, {x, y}]
SeedRandom[1]; NMinimize[rast, {x, y}, Method -> {"DifferentialEvolution", "PostProcess" -> False}]
SeedRandom[2]; NMinimize[rast, {x, y}, Method -> {"DifferentialEvolution", "PostProcess" -> False}]
NMinimize[rast, {x, y}, Method -> {"DifferentialEvolution", "PostProcess" -> False, "RandomSeed" -> 7}]
