---
source: src/pack.c
---
**Algorithm.** `FromPackedArray` is Mathematica's name (`Developer`​`FromPackedArray`) for
`FromNDArray` and is **the same C builtin** (`builtin_fromndarray`) registered under a second
name — an aliased builtin, not a DownValue that rewrites to the other, so it costs no extra
evaluation pass and cannot be shadowed. It undoes both forms of buffer storage: a packed
`List` becomes an ordinary `List` of separate element nodes (`NDArrayQ` then `False`), and a
visible `NDArray[...]` becomes the nested `List` of its entries (`ndarray_to_nested_list`).
Anything that is not an `EXPR_NDARRAY` is returned unchanged. It is the inverse of
`ToPackedArray`.

**Data structures.** Reads the dense row-major buffer and allocates one boxed `Expr` leaf per
element (`Integer`/`Real`/`True`/`False` by dtype), rebuilding the nested `List` from the
stored `int64` dims.

**Complexity / limits.** `O(n)`, allocating an expression node per element — the per-element
cost packing avoids, so this is where a large array stops being cheap to hold. The value is
unchanged; only the representation differs.
