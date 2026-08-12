# PseudoInverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PseudoInverse[m]`**

finds the Moore-Penrose pseudoinverse of a rectangular matrix m.

**`PseudoInverse[m, Tolerance -> t]`**

specifies that singular values smaller than t times the maximum singular value should be dropped.  With the default setting Tolerance -\> Automatic, the rationalisation precision of the input is used (Real -\> 53 bits, MPFR -\> input precision).

<details>
<summary>Notes</summary>

For non-singular square matrices m, the pseudoinverse coincides with the standard inverse: PseudoInverse\[m\] == Inverse\[m\]. PseudoInverse works on exact (Integer / Rational / Complex) matrices and on approximate (Real / MPFR) matrices.  For exact input the result is exact; for inexact input the input is rationalised, the pseudoinverse is computed in exact arithmetic via a full-rank decomposition, and the result is numericalised back to the input precision. Algorithm: row-reduce m to identify rank r and a full-rank decomposition m = B . C with B m x r and C r x n.  Then PseudoInverse\[m\] = ConjugateTranspose\[C\] . Inverse\[C . ConjugateTranspose\[C\]\] . Inverse\[ConjugateTranspose\[B\] . B\] . ConjugateTranspose\[B\]. When m is the zero matrix the pseudoinverse is the corresponding zero matrix of transposed shape.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

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

### Applications (5)

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

## Algorithm

matinv.c

Inverse and PseudoInverse.

```text
  Inverse[m]           -- exact / fraction-free Gauss-Jordan inversion
                          of a non-empty square matrix.  Lifted verbatim
                          from the previous src/linalg.c implementation.

  PseudoInverse[m]
  PseudoInverse[m,
      Tolerance -> t]  -- Moore-Penrose pseudoinverse of a rectangular
                          (or rank-deficient square) matrix.
```

Algorithm (PseudoInverse):

```text
  For an m x n matrix A with rank r > 0, compute the reduced row-echelon
  form R of A.  The first r non-zero rows of R, taken together as a
  r x n matrix C, span the row space of A; the columns of A at the pivot
  positions, taken together as a m x r matrix B, span the column space.
  This gives a full-rank decomposition A = B . C with rank(B) = rank(C) = r.
  The Moore-Penrose pseudoinverse is then

      A^+ = C^H . (C . C^H)^-1 . (B^H . B)^-1 . B^H

  where ^H is the conjugate transpose.  When A is invertible (m == n, r == n)
  the formula collapses to the standard inverse.  For the zero matrix
  (r == 0) the pseudoinverse is the n x m zero matrix.

  For inexact (Real / MPFR) matrices we rationalise the input at the
  minimum precision present (the common_rationalize_input pipeline used
  throughout the system), do every step in exact rational arithmetic so
  the rank is well-defined, then numericalise the final result back to
  that precision.  Tolerance -> Automatic uses the input precision.
```

Memory ownership follows the standard builtin contract: this file owns the `res` argument on success and frees it; on failure (returning NULL)

```text
the caller (evaluator) retains ownership.  Every intermediate matrix
```

allocated by Dot/Inverse/Transpose/Conjugate via eval_and_free is explicitly released.

## Implementation notes

**Algorithm.** `builtin_pseudoinverse` computes the Moore-Penrose pseudoinverse via an **exact full-rank decomposition**, not via SVD. `pseudoinverse_exact` row-reduces `A` (`mat_rref`), uses `find_pivots` to recover the rank `r` and the pivot columns, and forms a rank factorisation `A = B·C`: `B` is the `m × r` matrix of `A`'s pivot columns (`extract_columns`), `C` is the `r × n` matrix of the non-zero RREF rows (`extract_rows`). The pseudoinverse is then assembled from the closed form
`A⁺ = Cᴴ (C Cᴴ)⁻¹ (Bᴴ B)⁻¹ Bᴴ`,
using `hermitian_transpose`, `mat_mult`, and `mat_invert` on the small `r × r` Gram matrices, with the product finally `expand_matrix`-ed. The zero matrix (`r == 0`) returns the `n × m` zero matrix; an invertible square `A` collapses to the ordinary inverse.

**Data structures.** Standard `Expr*` `List`-of-`List` matrices throughout; the `r × r` Gram matrices `C Cᴴ` and `Bᴴ B` are full-rank by construction so `mat_invert` always succeeds. Inexact (`Real`/`MPFR`) input goes through the `common_rationalize_input` → exact-pipeline → `common_numericalize_result` round-trip at the input precision, giving a well-defined rank during row reduction.

**Limits.** A `Tolerance` option is parsed but currently has no effect on the exact pipeline. Non-rectangular input emits `PseudoInverse::matrix`.

- `Protected`.
- Works on rectangular and rank-deficient matrices over the integers,
  rationals, machine-precision reals, MPFR reals, exact complex
  (`Complex[a, b]` entries), and inexact complex.
- For non-singular square matrices, `PseudoInverse[m] == Inverse[m]`.
- Computes a full-rank decomposition `m = B . C` (with `B` `m x r` and
  `C` `r x n`) by row-reducing `m` to identify the rank `r` and the
  pivot columns, then returns
  `PseudoInverse[m] = ConjugateTranspose[C] . Inverse[C . ConjugateTranspose[C]] . Inverse[ConjugateTranspose[B] . B] . ConjugateTranspose[B]`.
- For the `m x n` zero matrix, returns the `n x m` zero matrix.
- **A machine-real matrix takes the SVD instead.** `A⁺ = V Σ⁺ Uᵀ` from one thin
  LAPACK `gesdd`, with the singular values at or below the cutoff zeroed. The
  rationalise-and-row-reduce pipeline described above is the right answer for an
  exact or symbolic matrix and an unusable one for a machine matrix — it did not
  finish on a 300 × 300 in 180 s, where the SVD path takes 23 ms. An exact
  matrix is unaffected and still answers exactly.
- MPFR matrices, and any `Real` matrix small enough not to be held as a
  [packed list](../packed-arrays/index.md), are rationalised at the input precision,
  computed exactly to preserve rank, and numericalised back.
- The `Tolerance` option accepts `Automatic` (default), a non-negative
  number, or a non-negative `Rational`. It takes effect on the SVD path, where
  there are singular values to truncate: `Automatic` is
  `max(m, n) × $MachineEpsilon × σ₁` — LAPACK's own rank criterion, and NumPy's
  default `rcond` — and an explicit `t` is a fraction of the largest singular
  value `σ₁`. The exact pipeline has no singular values and ignores it.
- Issues `PseudoInverse::matrix` warning and returns unevaluated if the
  argument is not a non-empty rank-2 tensor.
- Returns unevaluated when an unknown option name is supplied or
  `Tolerance` receives a negative value.
- Satisfies the Moore-Penrose identities
  `m . p . m == m` and `p . m . p == p` for `p = PseudoInverse[m]`.

**Attributes:** `Protected`.

## See also

[Rational](../../arithmetic/Rational/)

## References

- A. Ben-Israel and T. N. E. Greville, *Generalized Inverses: Theory and Applications*, 2nd ed. (Springer, 2003).
- Source: [`src/linalg/inv.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/inv.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_matinv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_matinv.c)
- Tests: [`tests/test_matlstsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_matlstsq.c)
- Tests: [`tests/test_ndarray_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_linalg.c)

## Notes & additional examples

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
