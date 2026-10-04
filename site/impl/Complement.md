---
source: src/list/setops.c
---
**Algorithm.** `builtin_complement` gives the sorted distinct elements of its
first argument that appear in *none* of the later operands (set difference),
using the head of the first argument (need not be `List`); all operands must
share that head. It is the structural twin of `Intersection` with the membership
test inverted — a candidate survives when it is absent from every later operand.
After locating a trailing `SameTest` option (`Automatic` means the default), the
default path is `O(total)`: deduplicate the first operand through the file-local
chained hash set, then drop any candidate that hits a later operand's hash set.
`SameTest -> f` switches to an `O(n^2)` path with `f[a, b] === True` as the
equivalence relation; unlike `Intersection` (which keeps the greatest member of a
class), `Complement` keeps the canonically-*smallest* member, so the first
operand is sorted ascending and the first element seen for each class is the
representative. The result is deduplicated and sorted with `expr_compare`.

**Data structures.** The same `HashTable` of chained `HashNode` buckets
(`expr_hash`/`expr_eq` keys) and an `Expr**` candidate array as `Intersection`;
the buffer fast path uses `int64_t*` arrays. Unlike `Union`/`Intersection`,
`Complement` is order-sensitive in its first argument, so it is *not* `Flat`/
`OneIdentity`.

**Complexity / limits.** `O(total)` default, `O(n^2)` with `SameTest`. A rank-1
buffer of exact integers (invisible packed `List` or explicit `NDArray[...]`)
takes `setop_packed`'s machine fast path — a sorted set-difference over `int64`
words (`acc \ b`), or direct range-indexing when bounded — so `Complement` is on
`pack.c`'s `AWARE` and `INT64_OK` lists and keeps its input's representation;
reals, a custom `SameTest`, or a non-`List` head fall to the general path. There
is **no** `Compile[]` lowering (`CompileDiagnostics` reports `Compiled -> False`).
`Protected`.
