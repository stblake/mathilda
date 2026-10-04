---
source: src/list/take_drop.c
---
**Algorithm.** `builtin_drop` is the complement of `Take`: it removes the
elements selected by one or more sequence specifications and keeps the rest. Each
spec is resolved by `get_seq_spec_indices` into a concrete 1-based index list —
supporting a plain count `n` (front) or `-n` (back), `{m, n}` (a contiguous
block), `{m, n, s}` (a strided slice), `{m}` (a single element), `All`, `None`,
and `UpTo[k]`. `apply_take_drop` then recurses: at each level it marks the
selected indices for removal and rebuilds the enclosing function from the
survivors, descending into successive specs for multi-dimensional drops (so
`Drop[matrix, {2}, {2}]` deletes a row and a column at once).

**Data structures.** A boolean `keep` mask per level plus the computed index
array; the result reuses the original head, deep-copying the surviving children.
An `NDArray` first argument takes the native buffer path `ndstruct_drop` — a
contiguous leading-axis slice on the machine buffer, no boxing.

**Complexity / limits.** `O(n)` per level in the length of that level. Indices are
1-based, negative indices count from the end, and an out-of-range or invalid spec
returns `NULL`, leaving the expression unevaluated.
