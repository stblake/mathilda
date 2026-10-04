---
source: src/list/setops.c
---
**Algorithm.** `builtin_intersection` gives the sorted list of elements common to
all operands (or, for a single argument, its sorted distinct elements), using the
head of the first argument (need not be `List`); every operand must share that
head. It first locates a trailing `SameTest -> f` option. With the default test
it is `O(total)`: deduplicate the first operand into a candidate array through the
file-local chained hash set (`ht_insert`/`ht_find`, keyed by `expr_hash`/
`expr_eq`), then for each later operand build a hash set of its members and keep
only candidates present in it. With `SameTest -> f` it switches to an `O(n^2)`
path that treats `f[a, b] === True` as the equivalence relation and keeps the
canonically-greatest member of each class. The surviving candidates are sorted
with `expr_compare` for the result.

**Data structures.** A `HashTable` of chained `HashNode` buckets for the default
path; a borrowed/`expr_copy`'d `Expr**` candidate array carried through the
intersection passes. The int64 buffer fast path uses plain `int64_t*` arrays.

**Complexity / limits.** `O(total)` with the default test, `O(n^2)` with
`SameTest`. A rank-1 buffer of exact integers (from either the invisible packed
`List` or an explicit `NDArray[...]`) takes the machine fast path
`setop_packed` — a sorted-merge intersection over `int64` words (or direct
range-indexing when the value range is bounded) — so `Intersection` is on
`pack.c`'s `AWARE` and `INT64_OK` lists and keeps whichever representation it was
given; `setop_any_nd` guards it, and reals, a custom `SameTest`, or a non-`List`
head take the general path. There is **no** `Compile[]` lowering
(`CompileDiagnostics` reports `Compiled -> False`). `Flat`, `OneIdentity`,
`Protected`.
