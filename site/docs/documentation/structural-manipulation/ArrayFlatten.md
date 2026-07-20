# ArrayFlatten

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ArrayFlatten[a]
    creates a single flattened matrix from a matrix of matrices, forming a
    block matrix. Blocks sharing a row must agree on their first dimension
    and blocks sharing a column on their second; elements shallower than a
    block (e.g. 0) are treated as scalars and replicated to fill.
ArrayFlatten[a, r]
    flattens out r pairs of levels of a rank-2r array, giving a rank-r array.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= m = {{1, 2}, {3, 4}}; ArrayFlatten[{{0, 0, m}, {m, m, 0}}]
Out[1]= {{0, 0, 0, 0, 1, 2}, {0, 0, 0, 0, 3, 4}, {1, 2, 1, 2, 0, 0}, {3, 4, 3, 4, 0, 0}}
```

## Implementation notes

- `Protected`. Default `r = 2`.
- Blocks must fit: matrices in the same grid row must share their first dimension, and matrices in the same column their second; the output size along an axis is the sum of the block sizes. Disagreeing blocks leave the call unevaluated.
- Elements whose array depth is less than `r` are treated as scalars and replicated to fill a rank-`r` block (e.g. `0` becomes a zero block).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
