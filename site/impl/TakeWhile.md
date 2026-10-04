---
source: src/funcprog.c
---
**Algorithm.** `builtin_takewhile` returns the leading run of elements satisfying
`crit`, keeping the collection's head. Like `LengthWhile` it first tries the
compiled-predicate fast path (`pred_run_length` finds the run length on the buffer,
then `pred_leading_slice` copies it, so neither the scan nor the result materialises
boxed nodes). A visible `NDArray` is materialised and repacked
(`ndstruct_delist_repack`). Otherwise `leading_run_length` locates the run and the
first `k` elements are copied into a fresh collection; for an association the values
are tested and the head is `Association`.

**Data structures.** A packed buffer on the compiled path; an `Expr**` copy of the
leading `k` elements on the interpreted path.

**Complexity / limits.** O(k) for a run of length `k` (the scan stops at the first
failure), plus O(k) to copy the slice. `LengthWhile` returns only the count.
