# NegativeDefiniteMatrixQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NegativeDefiniteMatrixQ[m]
    gives True if m is explicitly negative definite, and False
    otherwise.

A matrix m is negative definite if Re[Conjugate[x] . m . x] < 0
for every nonzero vector x.  Equivalently, -m is positive
definite, i.e. the negated Hermitian part has only positive
eigenvalues.  The test is performed by attempting a Cholesky
factorisation of -(m + ConjugateTranspose[m]) / 2; on numeric
matrices this is dispatched to BLAS/LAPACK's dpotrf (real) or
zpotrf (complex) when available.  Returns False on non-numeric,
non-square, ragged, empty, or higher-rank tensor inputs.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

**Algorithm.** `builtin_negative_definite_matrix_q` tests `Re[Conjugate[x].m.x] < 0` for all nonzero `x`, equivalently that `−m` is positive definite. It mirrors `PositiveDefiniteMatrixQ` exactly but negates entries at load time: after the square-matrix shape gate and the `(re, im)`-double coercion (`ndq_leaf_to_double`), it loads `−m` into a column-major buffer, forms the Hermitian part of `−m` in the upper triangle, checks its diagonal is strictly positive (i.e. `m`'s diagonal strictly negative), and runs Cholesky via LAPACK `dpotrf`/`zpotrf` (with an in-house fallback). `info == 0` on `−m`'s Hermitian part ⇔ `m` negative definite ⇒ `True`; non-numeric/non-coercible entries give `False`. Wrong arity emits `NegativeDefiniteMatrixQ::argx`.

- `Protected`.
- Equivalent to: `-m` is positive definite, i.e. the negated Hermitian

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/negdef_q.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/negdef_q.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
