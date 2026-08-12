# IdentityMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IdentityMatrix[n] gives the n x n identity matrix.`**

**`IdentityMatrix[{m, n}] gives the m x n identity matrix.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= IdentityMatrix[3]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}

In[2]:= IdentityMatrix[{2, 3}]
Out[2]= {{1, 0, 0}, {0, 1, 0}}
```

### Applications (3)

```mathematica
In[1]:= IdentityMatrix[3]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

A two-element argument gives a rectangular identity (1s on the main diagonal,
0s elsewhere):

```mathematica
In[1]:= IdentityMatrix[{2, 3}]
Out[1]= {{1, 0, 0}, {0, 1, 0}}
```

It is the multiplicative identity for matrix products — multiplying any matrix
by a conformant identity leaves it unchanged:

```mathematica
In[1]:= IdentityMatrix[4] . HilbertMatrix[4] == HilbertMatrix[4]
Out[1]= True
```

## Implementation notes

`builtin_identitymatrix` accepts either an integer `n` (square `n×n`) or a pair `{m, n}` of integers, and constructs a `List` of `List`s with `Integer` `1` on the main diagonal (`i == j`) and `0` elsewhere. Non-integer or malformed dimension specs are returned unevaluated (`expr_copy(res)`). The output is exact integer entries; no numeric or symbolic processing is involved.

- `Protected`.
- Generates exact integer outputs (`1` on main diagonal, `0` elsewhere).
- Will remain unevaluated if arguments are symbolic or negative.

**Attributes:** `Protected`.

## References

- Source: [`src/linalg/construct.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/construct.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_diagonal_matrix_q.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal_matrix_q.c)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)

## Notes & additional examples

### Notes

`IdentityMatrix[n]` gives the `n x n` identity; `IdentityMatrix[{m, n}]` gives the `m x n` rectangular identity. Use it as a seed for matrix algebra (e.g. `MatrixPower`, characteristic-matrix constructions) and as the neutral element of `Dot`.
