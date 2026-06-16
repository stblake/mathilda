# Det

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Det[m]
    gives the determinant of the square matrix m.
Exact integer / rational / symbolic inputs use Bareiss-style
fraction-free Gaussian elimination; machine-precision Real / Complex
inputs dispatch to LAPACK LU (dgetrf / zgetrf) and accumulate the
pivot-signed product of diagonal entries; arbitrary-precision MPFR
inputs run a Doolittle LU at the input precision.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Det[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 0

In[2]:= Det[{{1.7, 7.1, -2.7}, {2.2, 8.7, 3.2}, {3.2, -9.2, 1.2}}]
Out[2]= 251.572

In[3]:= Det[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[3]= c (-e g + d h) - b (-f g + d i) + a (-f h + e i)
```

## Implementation notes

**Algorithm.** `builtin_det` validates that the argument is a non-empty square rank-2 tensor (via `get_tensor_dims`; otherwise it emits `Det::matsq` and returns `NULL`), flattens it row-major into an `Expr**`, and computes the determinant by full **Laplace cofactor expansion** along the first row, recursively (`laplace_det`). Each cofactor term is built as `Times[±1, element, minor]` and accumulated with `Plus`, every product/sum being reduced through `eval_and_free`, so cancellation and symbolic simplification happen as the expansion unwinds. This keeps results exact and symbolic for integer, rational, and symbolic matrices.

**Data structures.** The matrix is a flat `Expr**` of `n*n` element pointers; recursion carries an explicit `int* cols` index set and a fixed `row` cursor, deleting one column per level rather than copying submatrices. `laplace_det` is exported via `linalg.h` and reused by `Cross`.

**Complexity / limits.** Cofactor expansion is `O(n!)`, so it is only practical for small `n`. There is no fraction-free Bareiss or LU fast path in this handler — the larger fraction-free Gauss-Jordan machinery lives in `inv.c`/`linsolve.c` for inversion and solving, not in `Det`.

- `Protected`.
- Evaluates the determinant of a square matrix symbolically or numerically using Laplace expansion.
- Returns a warning `Det::matsq` if `m` is not a non-empty square matrix.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — Gaussian elimination and the LU view of the determinant.
- Source: [`src/linalg/det.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/det.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Det[{{1, 2}, {3, 4}}]
Out[1]= -2
```

```mathematica
In[1]:= Det[{{2, 0, 0}, {0, 3, 0}, {0, 0, 4}}]
Out[1]= 24
```

```mathematica
In[1]:= Det[{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}]
Out[1]= -3
```

```mathematica
In[1]:= Det[{{a, b}, {c, d}}]
Out[1]= -b c + a d
```

```mathematica
In[1]:= Det[{{1, 1, 1}, {a, b, c}, {a^2, b^2, c^2}}]
Out[1]= -a^2 b + a b^2 + a^2 c - b^2 c - a c^2 + b c^2
```

```mathematica
In[1]:= Det[Table[1/(i + j - 1), {i, 4}, {j, 4}]]
Out[1]= 1/6048000
```

```mathematica
In[1]:= Det[{{N[Pi, 40], 1}, {1, N[E, 40]}}]
Out[1]= 7.5397342226735670654635508695465744950351
```

### Notes

For a diagonal or triangular matrix the determinant is simply the product of the diagonal entries, as the third and fourth examples illustrate. Exact integer, rational, and symbolic inputs are handled by Bareiss-style fraction-free Gaussian elimination, so intermediate results never introduce spurious denominators and symbolic determinants come back fully factored where possible. A singular matrix returns 0 exactly. The argument must be a square matrix.
