# Diagonal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Diagonal[m]`**

gives the list of elements on the leading diagonal of the matrix m (length Min\[rows, cols\], so it works for a non-square m).

**`Diagonal[m, k]`**

gives the elements on the k-th diagonal of m: k \> 0 above the leading diagonal, k \< 0 below.  An out-of-range k gives {}.

## Examples (11)

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

### Applications (6)

The leading diagonal

```mathematica
In[6]:= Diagonal[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[6]= {1, 5, 9}
```

The first superdiagonal, k > 0 above

```mathematica
In[7]:= Diagonal[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1]
Out[7]= {2, 6}
```

The first subdiagonal, k < 0 below

```mathematica
In[8]:= Diagonal[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, -1]
Out[8]= {4, 8}
```

Non-square: length is Min[rows, cols]

```mathematica
In[9]:= Diagonal[{{1, 2, 3}, {4, 5, 6}}]
Out[9]= {1, 5}
```

Symbolic entries are copied through verbatim

```mathematica
In[10]:= Diagonal[{{a, b}, {c, d}}]
Out[10]= {a, d}
```

A machine array rides the rank-2 buffer straight to rank-1

```mathematica
In[11]:= Diagonal[NDArray[{{1., 2.}, {3., 4.}}]]
Out[11]= NDArray[{1.0, 4.0}]
```

## Algorithm

Diagonal[m] / Diagonal[m, k] -- extract the k-th diagonal of a matrix.

Diagonal[m] gives the leading diagonal {m[[1,1]], m[[2,2]], ...} (length Min[rows, cols], so it works for a non-square m). Diagonal[m, k] gives the k-th diagonal: k > 0 above the leading diagonal, k < 0 below. An in-range but empty diagonal (|k| beyond the matrix) gives {}.

A machine-precision matrix (a packed List or a visible NDArray) takes the buffer fast path in ndstruct_diagonal (rank-2 buffer straight to rank-1). The generic path below walks the nested List, so it also handles any rank >= 2: the diagonal of a rank-n tensor is rank-(n-1), since each m[[i, i+k]] is itself an (n-1)-tensor and is copied through verbatim.

## Implementation notes

**Algorithm.** `builtin_diagonal` extracts the `k`-th diagonal of a matrix. `Diagonal[m]` gives the leading diagonal `{m[[1,1]], m[[2,2]], ...}` of length `Min[rows, cols]` (so it works for a non-square `m`); `Diagonal[m, k]` takes the `k`-th diagonal, `k > 0` above the leading one and `k < 0` below. From `k` it sets a `(start_row, start_col)` origin and walks `t = 0, 1, ...`, reading `m[[start_row + t, start_col + t]]` until a row index runs past the matrix, a column index runs past that row, or a row is not itself a `List`. A machine-precision argument (a packed `List` or a visible `NDArray`) is detected by `is_ndarray` and handed to `ndstruct_diagonal`; the generic path handles everything else and leaves a non-matrix (a vector, scalar, or non-list) unevaluated, as Mathematica does.

**Data structures.** The generic path walks the nested-`List` `Expr`, collecting `expr_copy`'d entries into a `cap`-doubling `Expr**` buffer and wrapping them in a `List`. Because each entry is copied verbatim, the diagonal of a rank-`n` tensor is a rank-`(n-1)` structure (each `m[[i, i+k]]` is itself an `(n-1)`-tensor). The buffer fast path `ndstruct_diagonal` (`src/ndstruct.c`) reads the rank-2 buffer directly: it `memcpy`s `esz` bytes per diagonal element straight into a fresh rank-1 buffer of the same `dtype` and returns `expr_new_ndarray_like(a, 1, ...)` — so a packed `List` stays packed, a visible `NDArray` stays visible, and the dtype is preserved with no boxing. `Diagonal` is on `src/pack.c`'s `AWARE` list, so the evaluator's transparency gate passes a packed argument through instead of materialising it.

**Complexity / limits.** `O(len)` element copies, where `len` is the diagonal length. An out-of-range `|k|` (beyond the matrix) gives `{}` — the generic loop is simply empty, and the buffer path returns early on `k <= -R || k >= C`. Non-square and higher-rank inputs are both supported; a non-matrix argument is left unevaluated.

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

## References

**See also:** [List](../../other-advanced/List/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/linalg/diagonal.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/diagonal.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_diagonal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal.c)
- Tests: [`tests/test_jordandecomp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_jordandecomp.c)

## Notes & additional examples

### Notes

`Diagonal[m, k]` indexes the diagonal by its offset `k` from the leading one:
`k > 0` runs above it, `k < 0` below. An `|k|` beyond the matrix gives `{}`.

Because each entry `m[[i, i+k]]` is copied unchanged, `Diagonal` works on more than
square numeric matrices: it handles non-square inputs (the length is `Min[rows, cols]`)
and higher-rank tensors (the diagonal of a rank-`n` tensor is a rank-`(n-1)` structure).

A packed `List` or a visible `NDArray` takes the buffer fast path, which `memcpy`s the
diagonal into a fresh array of the same dtype — a packed input stays packed, a visible
`NDArray` stays visible.
