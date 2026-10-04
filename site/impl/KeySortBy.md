---
source: src/assoc.c
---
**Algorithm.** `builtin_keysortby` computes `f[key]` once per entry, storing the
owned `f`-values in a parallel array, then runs a **stable insertion sort** keyed by
those values (compared with `expr_compare`). Stability is deliberate: entries whose
`f`-values tie keep their association order. The entries are rebuilt into a fresh
`Association` in the sorted order.

**Data structures.** Two parallel `Expr**` arrays — the entry copies and their
`f`-values — reordered together by the insertion sort.

**Complexity / limits.** O(n²) comparisons in the worst case (insertion sort, chosen
for a cheap stable order) with one `f` evaluation per entry, which is well suited to
the modest sizes associations usually reach. `KeySort` sorts by the keys themselves;
`SortBy` is the list/association analogue that sorts by a function of the values.
