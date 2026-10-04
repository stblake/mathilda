---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/assoc.c
---
**Algorithm.** `builtin_groupby` returns `<|f[x] -> {x, ...}|>`, preserving
first-key order. It evaluates the key function once per element, looks the key up
in a `KeyIndex`, and appends the element (a copy) to that key's growing group
buffer. A third argument reduces each group (`GroupBy[list, f, red]` →
`<|k -> red[group]|>`); a `keyfn -> valfn` second argument groups by `keyfn[x]`
but collects `valfn[x]`; and a `List` of key functions grooves multi-level
grouping into nested associations, the reducer applying only at the innermost
level. Over an `Association` the entries are grouped by `f[value]` with keys
preserved.

**Data structures.** A `KeyIndex` open-addressing hash set over the owned group
keys, parallel with per-group dynamic `Expr**` buffers that double on demand; the
groups become `List`s (or sub-associations) under the group keys.

**Complexity / limits.** `O(n)` evaluations of the key function plus `O(n)`
hashing and copying; multi-level grouping recurses through the evaluator (each
level a `GroupBy` call). Returns `NULL` unless the first argument is a `List` or
association. `f` is an arbitrary function, so there is no packed fast path.
