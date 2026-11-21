# Symbolic Integration -- ParallelMixedTower: Kummer split, conics, kernel substitutions
Integrate[1/(Sqrt[x] + Sqrt[x + 1]), x, Method -> "ParallelMixedTower"]
Integrate[1/Sqrt[x^2 + a x + b], x, Method -> "ParallelMixedTower"]
Integrate[Sqrt[1 + Log[x]]/(x Log[x]), x, Method -> "ParallelMixedTower"]
Integrate[Sqrt[1 + E^x], x, Method -> "ParallelMixedTower"]
Integrate[Sqrt[x + Sqrt[x^2 + 1]]/x, x, Method -> "ParallelMixedTower"]
