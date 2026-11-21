# Symbolic Integration -- Galois components and support groups
Integrate[Sqrt[x] Log[x] + Sqrt[x + 1] Log[x], x, Method -> "Bronstein"]
Check[Integrate[Sqrt[x] Log[x] + E^x/x, x, Method -> "Bronstein"], Last[$MessageList]]
