---
source: src/boolean.c
---
**Algorithm.** `Implies` is `Protected` (binary, no `Flat`/`Orderless`, since
material implication is neither associative nor commutative). `builtin_implies`
requires exactly two arguments and tests a short ladder of literal and structural
cases: `Implies[False, q]` and `Implies[p, True]` are `True`, `Implies[True, q]`
is `q`, `Implies[p, False]` is `Not[p]`, and `Implies[p, p]` (by `expr_eq`) is
`True`. Anything else returns `NULL` and stays symbolic.

**Data structures.** None; the two argument pointers are read in place and only
the chosen result is constructed.

**Complexity / limits.** `O(1)` — a fixed set of comparisons. No further
reasoning is done at this layer; `LogicalExpand` and `Reduce` rewrite
`Implies[p, q]` to its definition `!p || q` when a boolean normal form is wanted.
