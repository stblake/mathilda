# Symbolic Integration -- algebraic functions by Trager's algorithm
Integrate[1/(x^4 + 1)^(1/4), x]
Integrate[1/(x^3 + 1)^(1/3), x]
Integrate[(x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 1]), x]
Integrate[x/Sqrt[x^4 + 10 x^2 - 96 x - 71], x]
Integrate[(Sqrt[x] + (x + 1)^(1/3))/(x (x + 1)), x]
Integrate[Sqrt[x] (x + 1)^(1/3), x, Method -> "Trager"]
Check[Integrate[Sqrt[x] (x + 1)^(1/3), x, Method -> "Trager"], Last[$MessageList]]
