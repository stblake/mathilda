# LeanForm: the display translation into Lean 4 syntax.
LeanForm[(x + y)^2 == x^2 + 2 x y + y^2]
LeanForm[x/(y z) + 1/x - 3/4 a]
ToString[{Function[{x}, x + 1], Log[b], f[-3], Pi}, LeanForm]
ToString[a < b <= c && Not[p], LeanForm]
# The strict translator refuses whatever Lean would read differently.
Lean`Internal`ToLean[x/2, "Int", True, {}]
Lean`Internal`ToLean[1.5 x, "Rat", True, {}]
Lean`Internal`ToLean[Log[x], "Real", True, {}]
Lean`Internal`ToLean[x + 3/4, "Rat", True, {}, True]
