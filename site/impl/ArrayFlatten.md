---
source: src/list/array_flatten.c
---
**Algorithm.** `builtin_array_flatten` turns a rank-`r` *grid of blocks* into one
rank-`r` array — `ArrayFlatten[a]` is the `r = 2` block-matrix case, equal to
`Flatten[a, {{1, 3}, {2, 4}}]` (both spellings implemented directly here, as the
list-of-lists `Flatten` levelspec does not exist in this codebase).
`af_collect_grid` walks `a` down `r` `List` levels, recording the grid dimension
at each level and collecting a borrowed pointer to every block in row-major
order; a ragged grid or a non-`List` node makes it decline (`NULL`,
unevaluated). Each block is then classified by `af_tensor_dims`: a block whose
array depth is `>= r` is an *array block* whose first `r` dimensions give its
size along each grid axis, and anything shallower (an atom, a non-`List` head, a
too-shallow list) is a *scalar* replicated to fill its grid slot — this is how a
`0` becomes a zero block. Blocks must fit: all array blocks sharing a grid index
along an axis must agree on that dimension, or the call is left unevaluated. The
output extent along each axis is the sum of the fitted per-position lengths
(positions with no array block contribute `1`), and `af_build` lays out the
result by recursion over the output axes, mapping each output multi-index back to
its `(grid index, within-block index)` and copying the element (or the replicated
scalar).

**Data structures.** Fixed-size `int64_t` buffers of width `AF_MAX_R = 32` for
the grid dims `D`, output dims `O`, grid strides, and the scratch output index; a
growable `Expr**` of borrowed block pointers; per-axis `int64_t*` maps `map_i`
(output index → grid index) and `map_j` (output index → within-block index); and
an `int*` flagging each block array-vs-scalar. The result is freshly built nested
`List`s via `expr_copy` of the source elements.

**Complexity / limits.** `O(total output elements)` copies plus one classifying
pass over the blocks; `r` is capped at `AF_MAX_R`. A structural head over nested
`List`s — not on `pack.c`'s `AWARE` list and with no `Compile[]` lowering.
`Protected`.
