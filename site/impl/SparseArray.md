---
source: src/sparsearray.c
---
**Algorithm.** `SparseArray` is an **inert head** — Mathilda has no sparse storage.
There is no `builtin_sparsearray`: a `SparseArray[rules, dims, default, …]`
specification simply persists unevaluated, and `Head[SparseArray[{1 -> 5}, 3]]` is
`SparseArray`. All the real work lives in `sparse_array_to_dense`, which is reached
*only* through `Normal[SparseArray[…]]` (dispatched in `builtin_normal`,
`src/calculus/series.c`). It allocates one flat row-major buffer of `prod(dims)`
slots, fills it rule by rule (first writer wins), tops the untouched slots up with the
default (`0` unless given), and folds the buffer into nested `List`s. Accepted spec
forms: explicit position rules `{pos -> v}`, `{p1, p2} -> {v1, v2}`, pattern rules
such as `{i_, i_} -> 1` (which visit every slot and so need `dims`), `Band[start] ->
v`, a dense `List`, and the internal CSR form
`SparseArray[Automatic, dims, default, {1, {rowptr, colidx}, vals}]`.

**Data structures.** A transient dense `Dense{rank, dims[32], total, slots}` buffer of
`Expr*` (NULL = not yet written); no compressed/CSR representation is *kept* — the
object remains its symbolic specification. Positions are 1-based with negative
indices resolved from the end (`dense_offset`).

**Complexity / limits.** The `Normal` conversion is O(prod(dims)) in both time and
space, which is the honest cost: because there is no sparse backing store, arithmetic
and linear algebra are **not** accelerated on a `SparseArray` — it must be
`Normal`-ised to a dense array first. `Normal` refuses (leaving the call unevaluated)
above `SA_MAX_ELEMENTS = 2^27` slots or `SA_MAX_RANK = 32`, to avoid a multi-gigabyte
allocation. This is a deliberate, documented limitation.
