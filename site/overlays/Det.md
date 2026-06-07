---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — Gaussian elimination and the LU view of the determinant."
---
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
