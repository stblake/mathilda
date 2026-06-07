# Dot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
a . b . c or Dot[a, b, c]
    contracts the last index of each argument with the first index of
    the next: matrix-matrix, matrix-vector, vector-vector, and general
    tensor inner products.
Numeric machine-precision Real / Complex matrix-matrix dot dispatches
to BLAS dgemm / zgemm when available; exact and symbolic inputs use
the elementwise sum-of-products.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= {a, b, c} . {x, y, z}
Out[1]= a x + b y + c z

In[2]:= {{a, b}, {c, d}} . {x, y}
Out[2]= {a x + b y, c x + d y}

In[3]:= {{a, b}, {c, d}} . {{1, 2}, {3, 4}}
Out[3]= {{a + 3 b, 2 a + 4 b}, {c + 3 d, 2 c + 4 d}}
```

## Implementation notes

- `Flat`, `OneIdentity`, `Protected`.
- Contracts the last index in `a` with the first index in `b`.
- Applying `Dot` to a rank `n` tensor and a rank `m` tensor gives a rank `m+n-2` tensor.
- Scalar product of two vectors yields a scalar.
- Product of a matrix and a vector yields a vector.
- Product of two matrices yields a matrix.
- When arguments are not lists, `Dot` remains unevaluated.
- Gives an error message `Dot::dotsh` if the shapes of the inputs are not compatible.

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — matrix and vector products.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Dot[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]
Out[1]= {{19, 22}, {43, 50}}
```

```mathematica
In[1]:= {1, 2, 3} . {4, 5, 6}
Out[1]= 32
```

```mathematica
In[1]:= {{1, 2}, {3, 4}} . {5, 6}
Out[1]= {17, 39}
```

### Notes

`Dot` contracts the last index of each argument with the first index of the next, so it covers matrix-matrix products, matrix-vector products, and the vector dot product (which yields a scalar) with a single operator. The infix `.` form is equivalent to `Dot[...]` and chains for products of three or more arrays. Exact and symbolic inputs use the elementwise sum of products; machine-precision Real and Complex matrix-matrix products dispatch to a BLAS `gemm` kernel when available.
