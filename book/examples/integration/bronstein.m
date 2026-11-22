# Symbolic Integration -- a radical over a log/exp tower (Bronstein)
Integrate[(1 + 1/x)/(2 Sqrt[x + Log[x]]), x]
Integrate[Sqrt[x + Log[x]], x]
Check[Integrate[Sqrt[x + Log[x]], x], Last[$MessageList]]
Integrate[E^x/Sqrt[1 + E^x], x]
f = D[Sqrt[x - 1] (x + 1)^(1/3) E^x, x]
Integrate[f, x]
Integrate[Log[x]/Sqrt[1 + x^3], x]
Integrate[E^(x Sqrt[x^3 + 1])/x, x]
