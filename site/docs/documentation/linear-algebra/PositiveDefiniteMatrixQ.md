# PositiveDefiniteMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PositiveDefiniteMatrixQ[m]
    gives True if m is explicitly positive definite, and False
    otherwise.

A matrix m is positive definite if Re[Conjugate[x] . m . x] > 0
for every nonzero vector x.  Equivalently, the Hermitian part
(m + ConjugateTranspose[m]) / 2 has only positive eigenvalues.
The test is performed by attempting a Cholesky factorisation of
the Hermitian part; on numeric matrices this is dispatched to
BLAS/LAPACK's dpotrf (real) or zpotrf (complex) when available.
Returns False on non-numeric, non-square, ragged, empty, or
higher-rank tensor inputs.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{5, -1}, {-1, 4}}]
Out[1]= True

In[2]:= PositiveDefiniteMatrixQ[{{2.3, -1.2}, {0.6, 3.7}}]
Out[2]= True

In[3]:= PositiveDefiniteMatrixQ[{{1, 2 I}, {-I, 4}}]
Out[3]= True

In[4]:= PositiveDefiniteMatrixQ[{{Pi, -5, 2}, {E, -3, -3}, {5, Sqrt[2], 5}}]
Out[4]= False

In[5]:= PositiveDefiniteMatrixQ[{{1, a}, {b, 2}}]
Out[5]= False

In[6]:= PositiveDefiniteMatrixQ[Table[1/(i + j - 1), {i, 8}, {j, 8}]]
Out[6]= True
```

## Implementation notes

**Algorithm.** `builtin_positive_definite_matrix_q` tests whether `Re[Conjugate[x].m.x] > 0` for all nonzero `x`, equivalently that the Hermitian part `H = (m + ConjugateTranspose[m])/2` admits a Cholesky factorisation. After the square-matrix shape gate it coerces every entry to a `(re, im)` double (`pdq_leaf_to_double`: Integer, BigInt, Real, MPFR, exact `Rational`, `Complex`); any non-coercible leaf yields `False` (the predicate refuses claims it cannot prove). It builds `H` into the upper triangle of a column-major buffer (note: the input need not itself be Hermitian — only `H` matters), checks the diagonal is strictly positive, then runs Cholesky via LAPACK `dpotrf`/`zpotrf` (in-house `pdq_chol_real_inplace`/`pdq_chol_complex_inplace` fallback when `USE_LAPACK` is off). `info == 0` ⇔ positive definite ⇒ `True`. Wrong arity emits `PositiveDefiniteMatrixQ::argx`.

- `Protected`.
- Equivalent to: the Hermitian part `(m + ConjugateTranspose[m]) / 2`

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/posdef_q.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/posdef_q.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{2, 1}, {1, 2}}]
Out[1]= True
```

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{1, 2}, {2, 1}}]
Out[1]= False
```

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}]
Out[1]= True
```

```mathematica
In[1]:= PositiveDefiniteMatrixQ[{{2, I}, {-I, 2}}]
Out[1]= True
```

### Notes

`PositiveDefiniteMatrixQ[m]` tests whether the Hermitian part of `m` has only
positive eigenvalues, equivalently whether `Re[Conjugate[x] . m . x] > 0` for
every nonzero vector `x`. It is decided by attempting a Cholesky factorisation
of the Hermitian part. The second matrix `{{1,2},{2,1}}` has eigenvalues `3`
and `-1`, so it is indefinite and the test returns `False`. The third example
is the classic `3x3` second-difference (discrete Laplacian) matrix, which is
positive definite. The last example is a genuinely Hermitian complex matrix
`{{2, I}, {-I, 2}}` (note `Conjugate[I] = -I` off the diagonal) with
eigenvalues `1` and `3`; the complex path dispatches to LAPACK's `zpotrf` when
available. The function returns `False` on non-numeric, non-square, ragged, or
empty inputs.
