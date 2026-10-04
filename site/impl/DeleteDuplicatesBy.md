---
source: src/list/setops.c
---
**Algorithm.** `builtin_deleteduplicatesby` keeps the first element for each
distinct value of `f[element]`, preserving order. It walks the collection,
evaluates `f` on each element (on each *value*, for an association), and compares
the result against the `f`-values of the survivors so far; a new `f`-value adds
the element to the kept list, a repeat drops it. Over an association the surviving
entries are returned as an association with keys preserved.

**Data structures.** Two parallel arrays — the surviving elements/rules and their
owned `f`-values — with survivors compared directly by `expr_eq`. The survivor
count is typically small, so the linear scan of seen `f`-values stays within
budget; no hash index is built.

**Complexity / limits.** `O(n · k)` for `n` elements and `k` survivors
(`O(n²)` worst case when nearly everything is distinct), plus `n` evaluations of
`f`. Requires a 2-argument call with a non-atomic first argument; `f` is
arbitrary, so no packed fast path applies.
