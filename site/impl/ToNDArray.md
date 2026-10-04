---
source: src/pack.c
---
**Algorithm.** `builtin_tondarray` first strips an optional trailing `DataType -> "..."` rule
with `pack_take_dtype` (rightmost wins); after that there must be exactly one positional
argument. An argument that is already a packed `List` with no conflicting `DataType` is
returned by `expr_copy` (a no-op); a visible `NDArray[...]` is turned back into a nested
`List` first (`ndarray_to_nested_list`), so `ToNDArray` is also a way to say "same values, as
a `List`". The list is then packed by `pack_force_coerce` → `pack_build` with `min_elems = 0`
(no size threshold, unlike automatic packing) and `coerce = true`:

1. `pack_sniff` walks the nested list once, checking rectangularity and folding every leaf's
   class (`PK_INT`/`PK_REAL`/`PK_BOOL`) into one; with `coerce = true` a mix of machine
   `Integer` and `Real` folds to `PK_REAL` (the integers will widen to `double`), whereas the
   automatic path (`coerce = false`) declines that mix to keep `1 === 1.` observable.
2. The dtype is `int64` (all-integer), `float64` (real, or coerced mixed), or `bool`
   (all-`True`/`False`). An explicit `DataType` may *widen* an exact list to a float buffer
   but is refused when it would round a `Real` into an `int64` slot (`want == NDT_INT64 &&
   cls != PK_INT`), request a complex buffer, or mismatch bool-vs-number — in which case the
   original list is returned unchanged.
3. `pack_flatten` writes the leaves row-major into a freshly `malloc`'d buffer.

A list that is ragged, empty, or holds any `Rational`/`BigInt`/arbitrary-precision/`Complex`/
symbolic element simply comes back unchanged (never an error).

**Data structures.** The result is an `EXPR_NDARRAY` with `present_as = NDA_HEAD_LIST` — the
*packed-list* surface, so `Head` is still `List` and only `NDArrayQ`/`PackedArrayQ` reveal the
buffer. Storage is a dense row-major buffer of `ndt_elem_size(dt) * n` bytes carrying rank and
`int64` dims (up to `NDARRAY_MAX_RANK`).

**Complexity / limits.** `O(n)` in two passes (sniff then flatten) and one allocation. Ignores
the automatic-packing size threshold. Complex dtypes are rejected (no faithful round trip
yet); the packing contract is representation-only, so every downstream operation either
answers exactly as the ordinary list does or falls back to it.
