---
source: src/core.c
---
**Algorithm.** `builtin_level` collects the subexpressions of `e` that lie at the depths
selected by the level spec (an integer `n` = levels `1..n`, `Infinity`, `{n}` = exactly level
`n`, or `{min, max}`). It descends with the recursive `level_rec`, tracking the current depth
and appending matching nodes into a growable `Expr**` buffer, then wraps them in a `List`. The
option `Heads -> True` additionally includes function heads as level elements.
