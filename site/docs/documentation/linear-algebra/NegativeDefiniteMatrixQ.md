# NegativeDefiniteMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NegativeDefiniteMatrixQ[m]`**

gives True if m is explicitly negative definite, and False otherwise.

<details>
<summary>Notes</summary>

A matrix m is negative definite if Re\[Conjugate\[x\] . m . x\] \< 0 for every nonzero vector x.  Equivalently, -m is positive definite, i.e. the negated Hermitian part has only positive eigenvalues.  The test is performed by attempting a Cholesky factorisation of -(m + ConjugateTranspose\[m\]) / 2; on numeric matrices this is dispatched to BLAS/LAPACK's dpotrf (real) or zpotrf (complex) when available.  Returns False on non-numeric, non-square, ragged, empty, or higher-rank tensor inputs.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= NegativeDefiniteMatrixQ[{{-5, 1}, {1, -4}}]
Out[1]= True

In[2]:= NegativeDefiniteMatrixQ[{{-2.3, -1.2}, {0.6, -3.7}}]
Out[2]= True

In[3]:= NegativeDefiniteMatrixQ[{{-1, 2 I}, {-I, -4}}]
Out[3]= True

In[4]:= NegativeDefiniteMatrixQ[{{Pi, -5, 2}, {E, -3, -3}, {5, Sqrt[2], 5}}]
Out[4]= False

In[5]:= NegativeDefiniteMatrixQ[{{-1, a}, {b, -2}}]
Out[5]= False

In[6]:= NegativeDefiniteMatrixQ[Table[-1/(i + j - 1), {i, 8}, {j, 8}]]
Out[6]= True
```

### Applications (4)

```mathematica
In[7]:= NegativeDefiniteMatrixQ[{{-2, 0}, {0, 3}}]
Out[7]= False

In[8]:= NegativeDefiniteMatrixQ[{{-3, 1, 0}, {1, -3, 1}, {0, 1, -3}}]
Out[8]= True

In[9]:= NegativeDefiniteMatrixQ[{{-2, I}, {-I, -2}}]
Out[9]= True

In[10]:= NegativeDefiniteMatrixQ[-Table[1/(i + j - 1), {i, 3}, {j, 3}]]
Out[10]= True
```

## Algorithm

src/linalg/negdef_q.c

NegativeDefiniteMatrixQ -- the explicit negative-definiteness predicate.

A matrix m is negative definite iff Re[Conjugate[x] . m . x] < 0 for

```text
every nonzero vector x.  Equivalently, -m is positive definite, i.e.
```

the negated Hermitian part -(m + ConjugateTranspose[m]) / 2 has only positive eigenvalues and admits a Cholesky factorisation H = U^H U with a real positive diagonal.

For numeric input we load the matrix into a column-major double buffer, form the negated Hermitian part in-place into the upper triangle, and

```text
dispatch to LAPACK's `dpotrf` / `zpotrf`.  Cholesky returns info == 0
```

iff the operand is positive definite, so info == 0 here means the

```text
input was negative definite.  When USE_LAPACK is unavailable we fall
```

back to an in-house Cholesky.

For symbolic / non-numeric input we return False -- "explicitly negative definite" follows the Mathematica convention that the predicate refuses to make claims it cannot prove.

Diagnostics:

```text
  - argc != 1 -> `NegativeDefiniteMatrixQ::argx` to stderr, the call
    is left unevaluated (mirrors PositiveDefiniteMatrixQ / SquareMatrixQ).
```

Shape rejections that return False: non-list input, empty list, non-square (including ragged), and 3-D tensors.

## Implementation notes

**Algorithm.** `builtin_negative_definite_matrix_q` tests `Re[Conjugate[x].m.x] < 0` for all nonzero `x`, equivalently that `−m` is positive definite. It mirrors `PositiveDefiniteMatrixQ` exactly but negates entries at load time: after the square-matrix shape gate and the `(re, im)`-double coercion (`ndq_leaf_to_double`), it loads `−m` into a column-major buffer, forms the Hermitian part of `−m` in the upper triangle, checks its diagonal is strictly positive (i.e. `m`'s diagonal strictly negative), and runs Cholesky via LAPACK `dpotrf`/`zpotrf` (with an in-house fallback). `info == 0` on `−m`'s Hermitian part ⇔ `m` negative definite ⇒ `True`; non-numeric/non-coercible entries give `False`. Wrong arity emits `NegativeDefiniteMatrixQ::argx`.

- `Protected`.
- Equivalent to: `-m` is positive definite, i.e. the negated Hermitian
  part `-(m + ConjugateTranspose[m]) / 2` has only positive eigenvalues
  and admits a Cholesky factorisation with a real positive diagonal.
- On numeric matrices the test is performed by attempting Cholesky on
  the negated Hermitian part.  When `USE_LAPACK` is available the
  routine dispatches to BLAS/LAPACK's `dpotrf` (real) or `zpotrf`
  (complex); otherwise an in-house Cholesky is used.  Either returns
  `info == 0` iff the matrix is negative definite.
- Builds the Hermitian part regardless of whether the input is
  Hermitian, so e.g. a real matrix with non-symmetric entries is
  tested via `-(m + m^T) / 2`.
- For symbolic or otherwise non-coercible entries the predicate
  conservatively returns `False`; "explicitly negative definite" is
  not proved symbolically.
- Returns `False` (rather than leaving unevaluated) on non-matrix,
  non-square, ragged, empty, or higher-rank tensor inputs.
- Exactly one argument is accepted; any other count emits a
  Mathematica-compatible

  ```
  NegativeDefiniteMatrixQ::argx: NegativeDefiniteMatrixQ called with N arguments; 1 argument is expected.
  ```

  to `stderr` and leaves the call unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/linalg/negdef_q.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/negdef_q.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_ndarray_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_linalg.c)
- Tests: [`tests/test_negative_definite_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_negative_definite_matrix_q.c)

## Notes & additional examples

### Notes

`NegativeDefiniteMatrixQ[m]` tests whether `Re[Conjugate[x] . m . x] < 0` for
every nonzero `x`, equivalently that `-m` is positive definite. The diagonal
example fails because of the positive entry `3`. The symmetric tridiagonal matrix
is negative definite (eigenvalues all negative). The complex Hermitian matrix
`{{-2, I}, {-I, -2}}` shows that the test uses the Hermitian part. The last case
is the negated `3x3` Hilbert matrix — notoriously ill-conditioned yet still
detected as negative definite. The check is performed via an attempted Cholesky
factorisation of `-(m + ConjugateTranspose[m])/2`, dispatched to LAPACK on
numeric input.
