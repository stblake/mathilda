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

**Algorithm.** `builtin_dot` left-folds a chain `Dot[a, b, c, …]` by repeatedly contracting adjacent operand pairs through the `dot2` helper until no further contraction is possible. `dot2` contracts the last axis of `a` with the first axis of `b`: it reads both tensor shapes (`get_tensor_dims`), checks the shared dimension `K` matches (else `Dot::dotsh`), flattens both into row-major `Expr**` arrays, and for each output cell `(r, s)` forms `Plus[Times[a_rk, b_ks] for k]`, reducing each `Times`/`Plus` with `eval_and_free` so entries simplify symbolically. The result shape is `dimsA[0..rankA-2] ++ dimsB[1..rankB-1]`, rebuilt into nested `List`s by `build_tensor`. This single routine covers vector·vector (scalar), matrix·vector, matrix·matrix, and higher-rank tensor contractions uniformly. `dot2` is also reused by `MatrixPower` and the eigen solver.

**Data structures.** Dense flat `Expr**` buffers `flatA`/`flatB`/`flatC` in row-major order; dimension vectors are fixed `int64_t[64]`. If only one operand remains it is returned directly; if nothing contracted, `NULL` is returned (leave unevaluated). `Dot` is a thin chain driver — the actual numeric/symbolic arithmetic is delegated to the evaluator's `Plus`/`Times`.

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
- Source: [`src/linalg/dot.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/dot.c)
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
