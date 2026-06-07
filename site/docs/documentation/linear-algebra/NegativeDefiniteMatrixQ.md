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

- `Protected`.
- Equivalent to: `-m` is positive definite, i.e. the negated Hermitian

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
