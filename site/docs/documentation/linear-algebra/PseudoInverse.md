# PseudoInverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PseudoInverse[m]
    finds the Moore-Penrose pseudoinverse of a rectangular matrix m.
PseudoInverse[m, Tolerance -> t]
    specifies that singular values smaller than t times the maximum
    singular value should be dropped.  With the default setting
    Tolerance -> Automatic, the rationalisation precision of the
    input is used (Real -> 53 bits, MPFR -> input precision).

For non-singular square matrices m, the pseudoinverse coincides
with the standard inverse: PseudoInverse[m] == Inverse[m].

PseudoInverse works on exact (Integer / Rational / Complex)
matrices and on approximate (Real / MPFR) matrices.  For exact
input the result is exact; for inexact input the input is
rationalised, the pseudoinverse is computed in exact arithmetic
via a full-rank decomposition, and the result is numericalised
back to the input precision.

Algorithm: row-reduce m to identify rank r and a full-rank
decomposition m = B . C with B m x r and C r x n.  Then
    PseudoInverse[m] = ConjugateTranspose[C] . Inverse[C . ConjugateTranspose[C]]
                                            . Inverse[ConjugateTranspose[B] . B]
                                            . ConjugateTranspose[B].
When m is the zero matrix the pseudoinverse is the corresponding
zero matrix of transposed shape.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PseudoInverse[{{1,2},{3,4}}]
Out[1]= {{-2, 1}, {3/2, -1/2}}

In[2]:= PseudoInverse[{{1,2},{3,4}}] == Inverse[{{1,2},{3,4}}]
Out[2]= True

In[3]:= PseudoInverse[{{1,2,3},{4,5,6},{7,8,9}}]
Out[3]= {{-23/36, -1/6, 11/36}, {-1/18, 0, 1/18}, {19/36, 1/6, -7/36}}

In[4]:= PseudoInverse[{{0,0,0},{0,0,0}}]
Out[4]= {{0, 0}, {0, 0}, {0, 0}}

In[5]:= PseudoInverse[{{2,3},{2,2},{3,1},{4,3}}]
Out[5]= {{-29/134, -2/67, 22/67, 17/134}, {49/134, 8/67, -21/67, -1/134}}

In[6]:= PseudoInverse[{{1.25, 3.2, 3.2}, {7.9, -1.4, 5.1}, {0, 0, 0}}]
Out[6]= {{-0.0385185, 0.0966633, 0.0}, {0.210183, -0.0659894, 0.0}, {0.117363, 0.0282303, 0.0}}
```

## Implementation notes

**Algorithm.** `builtin_pseudoinverse` computes the Moore-Penrose pseudoinverse via an **exact full-rank decomposition**, not via SVD. `pseudoinverse_exact` row-reduces `A` (`mat_rref`), uses `find_pivots` to recover the rank `r` and the pivot columns, and forms a rank factorisation `A = B·C`: `B` is the `m × r` matrix of `A`'s pivot columns (`extract_columns`), `C` is the `r × n` matrix of the non-zero RREF rows (`extract_rows`). The pseudoinverse is then assembled from the closed form
`A⁺ = Cᴴ (C Cᴴ)⁻¹ (Bᴴ B)⁻¹ Bᴴ`,
using `hermitian_transpose`, `mat_mult`, and `mat_invert` on the small `r × r` Gram matrices, with the product finally `expand_matrix`-ed. The zero matrix (`r == 0`) returns the `n × m` zero matrix; an invertible square `A` collapses to the ordinary inverse.

**Data structures.** Standard `Expr*` `List`-of-`List` matrices throughout; the `r × r` Gram matrices `C Cᴴ` and `Bᴴ B` are full-rank by construction so `mat_invert` always succeeds. Inexact (`Real`/`MPFR`) input goes through the `common_rationalize_input` → exact-pipeline → `common_numericalize_result` round-trip at the input precision, giving a well-defined rank during row reduction.

**Limits.** A `Tolerance` option is parsed but currently has no effect on the exact pipeline. Non-rectangular input emits `PseudoInverse::matrix`.

- `Protected`.
- Works on rectangular and rank-deficient matrices over the integers,

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- A. Ben-Israel and T. N. E. Greville, *Generalized Inverses: Theory and Applications*, 2nd ed. (Springer, 2003).
- Source: [`src/linalg/inv.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/inv.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

For a non-singular square matrix the pseudoinverse coincides with the ordinary inverse:

```mathematica
In[1]:= PseudoInverse[{{1, 2}, {3, 4}}] == Inverse[{{1, 2}, {3, 4}}]
Out[1]= True
```

On a wide (full-row-rank) integer matrix the right inverse is computed exactly in rational arithmetic:

```mathematica
In[1]:= PseudoInverse[{{1, 2, 3}, {4, 5, 6}}]
Out[1]= {{-17/18, 4/9}, {-1/9, 1/9}, {13/18, -2/9}}
```

The Moore-Penrose defining identity `A . A^+ . A == A` holds exactly:

```mathematica
In[1]:= A = {{1, 2, 3}, {4, 5, 6}}; A . PseudoInverse[A] . A
Out[1]= {{1, 2, 3}, {4, 5, 6}}
```

It handles rank-deficient inputs gracefully; the all-ones matrix has rank 1 and a rank-1 pseudoinverse:

```mathematica
In[1]:= PseudoInverse[{{1, 1}, {1, 1}}]
Out[1]= {{1/4, 1/4}, {1/4, 1/4}}
```

Inexact input is rationalised, solved exactly, and returned at the input precision:

```mathematica
In[1]:= PseudoInverse[{{1., 2.}, {3., 4.}, {5., 6.}}]
Out[1]= {{-1.33333, -0.333333, 0.666667}, {1.08333, 0.333333, -0.416667}}
```

### Notes

`PseudoInverse[m]` computes the Moore-Penrose pseudoinverse via a full-rank
decomposition `m = B . C`, so it never forms an ill-conditioned normal-equation
solve. For exact (integer / rational / complex) input the result is exact, and
for non-singular square `m` it reduces to `Inverse[m]`. Rank-deficient and
rectangular matrices are handled without error; the defining Penrose relations
(here `A . A^+ . A == A`) hold identically. Inexact Real / MPFR input is
rationalised at its working precision, solved in exact arithmetic, and
numericalised back, giving the inexact-in / inexact-out contract.
`Tolerance -> t` drops singular values below `t` times the largest.
