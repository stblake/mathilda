---
source: src/core.c
---
**Algorithm.** `builtin_booleanq` is a one-argument predicate. It checks that its
single argument is an `EXPR_SYMBOL` whose interned name is either `SYM_True` or
`SYM_False`, and returns the symbol `True` or `False` accordingly. The argument is
evaluated first (ordinary, non-held evaluation), so `BooleanQ[2 > 1]` sees the already
reduced `True`. A wrong argument count routes through `builtin_arg_error`.

**Data structures.** Plain pointer comparison against the two interned boolean names;
no allocation beyond the returned result symbol.

**Complexity / limits.** `O(1)`. Attributes: `Protected`. Unlike `TrueQ` — which asks
"does this reduce to `True`?" and so maps everything that is not `True` to `False` —
`BooleanQ` tests the symbol itself, so `BooleanQ[False]` is `True` while `TrueQ[False]`
is `False`. Any non-boolean (a number, an unevaluated symbol, an unreduced
inequality) gives `False`.
