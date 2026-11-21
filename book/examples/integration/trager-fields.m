# Symbolic Integration -- Trager over Q(i), number fields, and Q(a)
Integrate[1/(x Sqrt[x^2 + I]), x]
Integrate[1/(x Sqrt[x^2 + Sqrt[2]]), x]
Integrate[1/((x - Sqrt[2]) Sqrt[x^2 + 1]), x]
Integrate[1/(x Sqrt[x^3 + a]), x]
Integrate[1/(x Sqrt[x^3 + a x + 1]), x, Method -> "Trager"]
Integrate[1/Sqrt[x^3 + a], x, Method -> "Trager"]
Integrate[(3 x^2 + a)/(2 Sqrt[x^3 + a x + b]), x]
