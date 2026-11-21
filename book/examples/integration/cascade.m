# Symbolic Integration -- the cascade, its budgets, and what is still open
{$TragerCascadeBudget, $BronsteinCascadeBudget}
Check[Integrate[Log[x]/Sqrt[1 + x^3], x], Last[$MessageList]]
Check[Integrate[Log[1 + Sqrt[Log[x]]]/x, x, Method -> "Bronstein"], Last[$MessageList]]
Integrate[Log[1 + Sqrt[Log[x]]]/x, x]
Check[Integrate[Sqrt[Log[x] + Sqrt[Log[x]]], x, Method -> "Bronstein"], Last[$MessageList]]
FreeQ[Integrate[Sqrt[Tan[x]], x, Method -> "Bronstein"], Root]
