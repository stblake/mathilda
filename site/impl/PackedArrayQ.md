---
source: src/ndarray.c
---
**Algorithm.** `builtin_packedarrayq` is a one-argument predicate that returns the symbol
`True` when its argument is a packed array — a `List` the system stores as a dense
machine-precision buffer — and `False` otherwise. The test is exactly `is_packed_list`, so it
is deliberately **narrower** than `NDArrayQ`: `NDArrayQ` is `True` for either surface of an
`EXPR_NDARRAY` (the packed-`List` surface *or* a visible `NDArray[...]` object), whereas
`PackedArrayQ` is `True` only for the packed-`List` surface. A visible `NDArray[...]` is a
distinct atom (`AtomQ` true, `ListQ` false) rather than a `List` that happens to be packed, so
it gives `False` — matching Mathematica (`Developer`​`PackedArrayQ`), which has no visible
`NDArray` head at all. The two predicates differ only on that one input.

For the answer to be correct the head must itself be **packed-aware**: `PackedArrayQ` is on
both the `AWARE` and `INT64_OK` lists in `src/pack.c`. Otherwise the evaluator's transparency
gate would materialise a packed argument into a plain `List` before the builtin runs, and it
could only ever observe an already-unpacked list and answer `False` for every packed input.

**Data structures.** None of its own — it inspects the argument's tag (`EXPR_NDARRAY` with the
packed-`List` presentation) and allocates only the returned `True`/`False` symbol.

**Complexity / limits.** `O(1)`. `Protected`.
