# ArrayFlatten

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArrayFlatten[a]`**

creates a single flattened matrix from a matrix of matrices, forming a block matrix. Blocks sharing a row must agree on their first dimension and blocks sharing a column on their second; elements shallower than a block (e.g. 0) are treated as scalars and replicated to fill.

**`ArrayFlatten[a, r]`**

flattens out r pairs of levels of a rank-2r array, giving a rank-r array.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= m = {{1, 2}, {3, 4}}; ArrayFlatten[{{0, 0, m}, {m, m, 0}}]
Out[1]= {{0, 0, 0, 0, 1, 2}, {0, 0, 0, 0, 3, 4}, {1, 2, 1, 2, 0, 0}, {3, 4, 3, 4, 0, 0}}
```

### Applications (3)

A 0 scalar fills a zero block

```mathematica
In[2]:= m = {{1, 2}, {3, 4}}; ArrayFlatten[{{0, 0, m}, {m, m, 0}}]
Out[2]= {{0, 0, 0, 0, 1, 2}, {0, 0, 0, 0, 3, 4}, {1, 2, 1, 2, 0, 0}, {3, 4, 3, 4, 0, 0}}
```

Two blocks sharing a grid row glue side by side

```mathematica
In[3]:= ArrayFlatten[{{IdentityMatrix[2], {{5}, {6}}}}]
Out[3]= {{1, 0, 5}, {0, 1, 6}}
```

A 2x2 grid of 2x2 blocks becomes a 4x4 matrix

```mathematica
In[4]:= ArrayFlatten[{{m, m}, {m, m}}]
Out[4]= {{1, 2, 1, 2}, {3, 4, 3, 4}, {1, 2, 1, 2}, {3, 4, 3, 4}}
```

## Implementation notes

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

- `Protected`. Default `r = 2`.
- Blocks must fit: matrices in the same grid row must share their first dimension, and matrices in the same column their second; the output size along an axis is the sum of the block sizes. Disagreeing blocks leave the call unevaluated.
- Elements whose array depth is less than `r` are treated as scalars and replicated to fill a rank-`r` block (e.g. `0` becomes a zero block).

**Attributes:** `Protected`.

## References

- Source: [`src/list/array_flatten.c`](https://github.com/stblake/mathilda/blob/main/src/list/array_flatten.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_array_flatten.c`](https://github.com/stblake/mathilda/blob/main/tests/test_array_flatten.c)

## Notes & additional examples

### Notes

`ArrayFlatten[a]` treats `a` as a rank-2 grid of matrix blocks and glues them
into one matrix, as `MatrixForm[a]` would show them — equivalent to
`Flatten[a, {{1, 3}, {2, 4}}]`. Blocks must fit: matrices in the same grid row
share their first dimension and those in the same column share their second, and
the output size along an axis is the sum of the block sizes; disagreeing blocks
leave the call unevaluated. A block shallower than a matrix (an atom such as `0`)
is a scalar replicated to fill the rank-2 block its position demands — that is
how a `0` becomes a zero block. `ArrayFlatten[a, r]` flattens `r` level pairs of
a rank-`2r` array into a rank-`r` array (default `r = 2`).
