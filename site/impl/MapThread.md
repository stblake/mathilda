---
source: src/funcprog.c
---
**Algorithm.** `builtin_mapthread` applies `f` across *k* lists in parallel:
`MapThread[f, {{a1,a2,...}, {b1,b2,...}, ...}]` gives `{f[a1,b1,...],
f[a2,b2,...], ...}`. The outer container must be a `List` of the *k* expressions
to thread, and `f` takes *k* arguments. At the default level 1 the
`nd_mapthread2` fast path handles the common two-list case over packed buffers;
otherwise `mapthread_rec` descends `level` list levels in lock-step, building one
`f[...]` leaf per tuple of corresponding parts, then the whole result is
`evaluate`d once so `f`'s attributes fire and numeric leaves reduce. A structural
mismatch (the lists are not the same shape down through `level`) makes the
recursion return `NULL`, leaving the call unevaluated.

Any `NDArray` entry is materialised to a nested list before threading, so an
array threads exactly like the corresponding list; if an input was packed, the
result is repacked with that input's dtype (`map_try_repack`) — packed in, packed
out. `MapThread[f, {}]` is `{}`.

**Data structures.** The *k* thread entries are borrowed from the argument list
(copied into a scratch array only when an NDArray entry must be delisted).
`mapthread_rec` works on the borrowed sub-expression arrays, so no spine is copied
beyond the result being built.

**Complexity / limits.** `O(k · N)` leaf constructions for *N* tuples, plus one
final evaluation pass. The level argument must be a non-negative integer; the
entries must agree in shape down through that level or the call declines.
