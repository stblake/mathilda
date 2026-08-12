# Transpose

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Transpose[list]`**

Transposes the first two levels of list (swaps rows and columns of a matrix).

**`Transpose[list, {n1, n2, ...}]`**

Gives the transpose of list so that level k in list is level nk in the result. The spec must be a permutation of {1, ..., r} where r is the depth of list. A repeated index (e.g. {1, 1}) selects the corresponding diagonal. list must be a rectangular array.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Transpose[{{a, b}, {c, d}}]
Out[1]= {{a, c}, {b, d}}

In[2]:= Transpose[{{a, b}, {c, d}}, {1, 1}]
Out[2]= {a, d}
```

### Applications (5)

```mathematica
In[1]:= Transpose[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= {{1, 4}, {2, 5}, {3, 6}}
```

```mathematica
In[1]:= Transpose[{{1, 2}, {3, 4}}]
Out[1]= {{1, 3}, {2, 4}}
```

```mathematica
In[1]:= Transpose[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1}]
Out[1]= {1, 5, 9}
```

The Gram matrix `A^T . A` of a symbolic `2x3` design matrix is symmetric by
construction:

```mathematica
In[1]:= Transpose[{{a, b, c}, {d, e, f}}] . {{a, b, c}, {d, e, f}}
Out[1]= {{a^2 + d^2, a b + d e, a c + d f}, {a b + d e, b^2 + e^2, b c + e f}, {a c + d f, b c + e f, c^2 + f^2}}
```

For an antisymmetric matrix `M = -M^T`, the sum `M + Transpose[M]` vanishes:

```mathematica
In[1]:= m = {{0, 1, 2}, {-1, 0, 3}, {-2, -3, 0}}; m + Transpose[m]
Out[1]= {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Transpose then Dot (fused?) | 37.1 s | 63.1 s | 18.6 s |
| Partition window 8, offset 1 | 2.39 s | 32.7 s | 5.24 s |
| Transpose 2000x2000 | 1.58 s | 0.736 s | 2.73 s |
| Take rows 1;;1000 of 2000x2000 | 0.215 s | 0.275 s | 0.213 s |
| column slice m[[All, 1]] | 0.004 s | 0.007 s | 0.001 s |
| ArrayReshape 2x10^6 to 1000x2000 | -- | 0.119 s | 0.216 s |

## Implementation notes

**Algorithm.** `builtin_transpose` swaps the levels of a rectangular nested-`List` array. It
measures the array shape with `get_array_dimensions` (requiring depth ≥ 2 and rectangularity),
then either uses the default permutation `{2, 1, 3, …}` (one-argument form swaps the first two
levels) or the explicit permutation given as the second argument. `build_transposed` recursively
materialises the output array by mapping each output index path back to an input index path
through the permutation and copying the leaf via `get_element_at`. For a 2-D matrix (list of
rows) this is the ordinary `m[i][j] -> m[j][i]` swap. Returns `NULL` (unevaluated) for
non-rectangular or non-`List` inputs. `ConjugateTranspose` is `Conjugate[Transpose[...]]`.

- `Protected`.
- Works only on rectangular arrays.
- `Transpose[m, {1, 1}]` extracts the diagonal of a square matrix.

**Attributes:** `Protected`.

## References

- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — the matrix transpose and index permutations of tensors.
- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_conjugate_transpose.c`](https://github.com/stblake/mathilda/blob/main/tests/test_conjugate_transpose.c)
- Tests: [`tests/test_fit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fit.c)
- Tests: [`tests/test_hankelmatrix.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hankelmatrix.c)

## Notes & additional examples

### Notes

With one argument `Transpose` swaps the first two levels of a list, turning a `2x3` matrix into a `3x2` one. The optional permutation spec generalises this to arbitrary index reorderings of a rectangular array. A repeated index in the spec — `{1, 1}` in the third example — extracts the corresponding diagonal, here the main diagonal `{1, 5, 9}` of the `3x3` matrix. The spec must be a permutation of `{1, ..., r}` where `r` is the depth of the list.
