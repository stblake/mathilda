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

- `Protected`.
- Evaluates the determinant of a square matrix symbolically or numerically using Laplace expansion.
- Returns a warning `Det::matsq` if `m` is not a non-empty square matrix.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — Gaussian elimination and the LU view of the determinant.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
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

### Notes

For a diagonal or triangular matrix the determinant is simply the product of the diagonal entries, as the third and fourth examples illustrate. Exact integer, rational, and symbolic inputs are handled by Bareiss-style fraction-free Gaussian elimination, so intermediate results never introduce spurious denominators and symbolic determinants come back fully factored where possible. A singular matrix returns 0 exactly. The argument must be a square matrix.
