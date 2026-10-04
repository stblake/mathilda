---
source: src/pack.c
---
**Algorithm.** `builtin_fromndarray` undoes buffer storage. It takes one argument: anything
that is not an `EXPR_NDARRAY` (as reported by `is_ndarray`, which is true for *both* the
packed-`List` surface and a visible `NDArray[...]`) is returned unchanged by `expr_copy`;
otherwise the buffer is expanded into a nested `List` of separate element nodes by
`ndarray_to_nested_list`. So a packed `List` becomes an ordinary `List` (`NDArrayQ` then
`False`), and an `NDArray[...]` becomes the nested `List` of its entries. It is the inverse of
`ToNDArray` and does the same job as `Normal` on these two forms.

**Data structures.** Reads the dense row-major buffer and rebuilds one boxed `Expr` leaf per
element (`Integer` from an `int64` buffer, `Real` from `float64`, `True`/`False` from `bool`),
reassembling the nested `List` from the stored `int64` dims.

**Complexity / limits.** `O(n)` — it allocates an expression node per element, which is
exactly the per-element cost that packing exists to avoid, so this is the point at which a
large array stops being cheap to hold. Pure representation change: the values are identical to
the buffer's.
