# Diagonal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Diagonal[m]`**

gives the list of elements on the leading diagonal of the matrix m (length Min\[rows, cols\], so it works for a non-square m).

**`Diagonal[m, k]`**

gives the elements on the k-th diagonal of m: k \> 0 above the leading diagonal, k \< 0 below.  An out-of-range k gives {}.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[1]= {a, e, i}

In[2]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, 1]
Out[2]= {b, f}

In[3]:= Diagonal[{{a, b, c}, {d, e, f}, {g, h, i}}, -1]
Out[3]= {d, h}

In[4]:= Diagonal[{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}]
Out[4]= {1, 6, 11}

In[5]:= Diagonal[{{a, b, c, d}, {e, f, g, h}, {i, j, k, l}, {m, n, o, p}}, -2]
Out[5]= {i, n}
```

## Algorithm

Diagonal[m] / Diagonal[m, k] -- extract the k-th diagonal of a matrix.

Diagonal[m] gives the leading diagonal {m[[1,1]], m[[2,2]], ...} (length Min[rows, cols], so it works for a non-square m). Diagonal[m, k] gives the k-th diagonal: k > 0 above the leading diagonal, k < 0 below. An in-range but empty diagonal (|k| beyond the matrix) gives {}.

A machine-precision matrix (a packed List or a visible NDArray) takes the buffer fast path in ndstruct_diagonal (rank-2 buffer straight to rank-1). The generic path below walks the nested List, so it also handles any rank >= 2: the diagonal of a rank-n tensor is rank-(n-1), since each m[[i, i+k]] is itself an (n-1)-tensor and is copied through verbatim.

## Implementation notes

- `Protected`.
- `Diagonal[m]` gives the leading diagonal `{m[[1,1]], m[[2,2]], ...}`, of length
  `Min[rows, cols]` — so it works for a non-square `m`.
- `Diagonal[m, k]` gives the `k`-th diagonal: `k > 0` above the leading diagonal
  (superdiagonals), `k < 0` below (subdiagonals). An out-of-range `k` gives `{}`.
- Elements are copied verbatim, so symbolic, exact, machine and complex entries
  all flow through and an exact-integer diagonal stays exact.
- Generalises to higher rank: the diagonal of a rank-`n` array is rank-`(n-1)`,
  since each `m[[i, i+k]]` is itself an `(n-1)`-tensor.
- Machine-precision matrices (a packed `List` or a visible `NDArray`) slice the
  diagonal directly off the flat buffer (`ndstruct_diagonal`, `src/ndstruct.c`),
  and the rank-2 case lowers inside `Compile[]` (`ND_FNS` / `A_NDFN`).

**Attributes:** `Protected`.

## See also

[List](../../other-advanced/List/), [NDArray](../../linear-algebra/NDArray/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_diagonal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal.c)
