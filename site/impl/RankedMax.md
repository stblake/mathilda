---
source: src/sort.c
---
**Algorithm.** `builtin_ranked_max` selects the `n`-th largest element of a list
— `RankedMax[list, n]`, with `RankedMax[list, -n]` the `n`-th smallest, so
`RankedMax[list, n]` is exactly `RankedMin[list, -n]` (`RankedMax[list, 1]` is
`Max[list]`, `RankedMax[list, -1]` is `Min[list]`). It validates an integer,
nonzero `n` and delegates to the shared `ranked_select` core with `is_max =
true`, which negates `n` before computing the 1-based ascending rank `r`. As for
`RankedMin`, one scan picks the exact `expr_compare` path when every element is a
real numeric atom, else the `double`-key path (`±HUGE_VAL` for `±Infinity`, the
value of a real atom, else the machine-precision numericalisation of a symbolic
real); a non-real element makes it decline. `ranked_select_idx` quickselects the
`r`-th position (`O(m)` average, Hoare partition, stable index tiebreak) and the
exact element there is returned.

**Data structures.** The same `ranked_select` machinery as `RankedMin`: a
`size_t* idx` quickselected in place, an optional `double* keys`, and a
`RankedCtx` read by `ranked_cmp`; the result is an `expr_copy` of the selected
element in its exact form.

**Complexity / limits.** `O(m)` average quickselect. An `NDArray` argument takes
the buffer order-statistic path `ndred_ranked_max` (`int64` exactly, reals via an
`O(m)` quickselect) — `RankedMax` is on `pack.c`'s `AWARE` and `INT64_OK` lists.
It is also compilable:
`CompileDiagnostics[{{v, _Real, 1}}, RankedMax[v, 2]]` reports `Compiled -> True`.
`Protected`.
