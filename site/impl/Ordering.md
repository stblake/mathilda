---
source: src/sort.c
---
**Algorithm.** `builtin_ordering` gives the permutation of 1-based positions that
sorts `list`, so that `list[[Ordering[list]]] === Sort[list]`. It collects
borrowed subject pointers — for an `Association` each entry's *value* (the Rule's
second argument), otherwise the element itself — and sorts an index array `idx`.
Without a custom ordering function it runs libc `qsort` through
`ordering_index_compare`, which orders by `expr_compare` of the pointed-at
subjects and breaks ties by the original index ascending; that tiebreak is what
turns the unstable `qsort` into a **stable** argsort, so equal elements keep
their input order (`Ordering[{2, 2, 1}]` is `{3, 1, 2}`, and `Ordering[list, 1]`
names the first minimum). A custom `p` (3-arg form) uses the merge-sort
permutation `p_sort_perm`, placing ties exactly as `Sort[list, p]` does.
`Ordering[list, seq]` emits only a `Take`-style slice of the permutation
(`get_seq_spec_indices` handles an `Integer n`/`-n`, a `{m, n[, s]}` span,
`UpTo[k]`, or `All`). The result is always a `List` of `Integer` positions,
regardless of `list`'s head.

**Data structures.** A borrowed `Expr** subjects`, an `int64_t* idx` permutation,
and an optional `int64_t* sel` of selected ranks; the comparator reads a
file-static `(ordering_subjects, ordering_p)` context saved/restored around each
sort so a nested `Ordering` recurses correctly. The result is offered to the
packer (`pack_offer`), since a list of positions packs and is often used as
`Part` indices next.

**Complexity / limits.** `O(n log n)` comparisons. A packed or visible rank-1
`NDArray` takes the buffer fast path `ndstruct_ordering`, which argsorts the
machine words directly and returns a packed `int64` permutation — so, unlike
`Sort`, its result dtype is always `Integer` (it is on `pack.c`'s `AWARE` and
`INT64_OK` lists); rank ≥ 2, a complex dtype, or a custom comparator decline by
materialising and re-evaluating. It is also compilable:
`CompileDiagnostics[{{v, _Real, 1}}, Ordering[v]]` reports `Compiled -> True`
with `ResultType -> Array`, lowering to a delegated buffer argsort.
`Protected`.
