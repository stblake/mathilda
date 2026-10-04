---
source: src/list/nearest.c
---
**Algorithm.** `builtin_nearest` (two-argument form) returns the element(s) of a
list closest to a target `x`. `nearest_distance` forms `Abs[element - x]` by
composing `internal_subtract` and `internal_abs` through `eval_and_free`, so a
complex element uses its modulus for free. The shape follows `MinimalBy`, not a
quickselect: one pass finds the minimum distance, a second collects *every*
element whose distance equals it, in input order — so ties are returned together
(`Nearest[{1, 5, 10}, 3]` is `{1, 5}`) and the tie handling falls out of the
ascending collect pass rather than being coded.

**Numeric gate.** Every distance must be a real number (checked by
`list_real_number_q`), or the whole call declines — this covers a symbolic
element, a symbolic target, and a non-real complex in one test, and even rejects
a symbolic real such as `Pi`. The minimum and ties are decided with
`list_numeric_cmp`; an undecidable comparison declines rather than guessing.
`Nearest[{}, x]` is `{}` (checked before the gate).

**Complexity / limits.** O(n) distance evaluations and O(n) comparisons, with
O(n) peak extra memory; two evaluate passes per element dominate, so this is an
interpreter-speed path, not a buffer one — `Nearest` is deliberately *not* on
`pack.c`'s `AWARE` list, so a packed argument is materialised and a visible
`NDArray` (not a `List`) is left unevaluated rather than silently truncated. Only
the two-argument form is implemented (no n-nearest, radius, rule or
`DistanceFunction` forms). `ATTR_PROTECTED`.
