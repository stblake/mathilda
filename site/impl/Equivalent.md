---
source: src/boolean.c
---
**Algorithm.** `Equivalent` is `Flat, Orderless, OneIdentity, Protected`, so
nesting is flattened and arguments are sorted before `builtin_equivalent` runs.
`Equivalent[]` and `Equivalent[e]` are trivially `True`. Otherwise one pass
records whether a literal `True` and/or `False` is present and collects the
distinct non-literal atoms (duplicates dropped by `expr_eq`). A `True` together
with a `False` forces `False`. If a single truth value is present it forces every
remaining atom — `True` leaves each atom as-is, `False` wraps each in `Not` — and
the results are joined with `And` (a lone atom is returned directly). With no
literals but duplicates removed, a smaller `Equivalent` is rebuilt; distinct
symbolic atoms with nothing to simplify return `NULL` (stay symbolic).

**Data structures.** A single `malloc`'d array of borrowed pointers to the
distinct atoms, plus two `bool`s for the literal flags; only the surviving atoms
are deep-copied into the `And`/`Equivalent`/`Not` result.

**Complexity / limits.** The dedup check is pairwise (`O(n²)` `expr_eq`), adequate
for small boolean expressions. Simplification is structural; the full
truth-functional meaning is left to `LogicalExpand`/`Reduce`/`FindInstance`, which
expand `Equivalent[a1, …, an]` to the cyclic conjunction
`Implies[a1, a2] && … && Implies[an, a1]`.
