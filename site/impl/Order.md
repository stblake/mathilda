---
source: src/sort.c
---
**Algorithm.** `builtin_order` is the user-facing surface of `expr_compare`, the
canonical comparator every sorting routine (`Sort`, `SortBy`, `OrderedQ`,
`Ordering`, ...) is built on. `Order[e1, e2]` requires exactly two arguments and
returns `1` if `e1` is before `e2` in canonical order, `-1` if after, `0` if
identical. `expr_compare` returns negative when `e1` sorts first, so the sign is
simply inverted. Comparison is **structural**, not by numerical value:
`Order[6, Pi]` is `1` (the `Integer` `6` sorts before the symbol `Pi`) whereas
`Order[6, N[Pi]]` is `-1` (two numeric atoms, compared by value) — the canonical
order ranks reals by value, then strings, then symbols, then expressions by
length/head/parts.

**Data structures.** None of its own — one `expr_compare` call over the two
borrowed argument trees, returning an `EXPR_INTEGER` in `{1, 0, -1}`.

**Complexity / limits.** `O(min size)` of the two trees in the worst case (the
comparator short-circuits on the first difference). It **is** compilable:
`CompileDiagnostics[{{x, _Real}}, Order[x, 1]]` reports `Compiled -> True` with
`ResultType -> Integer`, lowering over machine numbers to `Sign[e2 - e1]`
(matching the interpreter's integer head); complex or array arguments fall back
to the interpreter. No packed/NDArray path — it is a scalar two-argument
comparator, not an element-wise map. `Protected`.
