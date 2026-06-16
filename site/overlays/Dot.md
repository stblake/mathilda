---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — matrix and vector products."
---
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

```mathematica
In[1]:= {{1, 2}, {3, 4}} . {{0, 1}, {1, 0}} . {{1, 2}, {3, 4}}
Out[1]= {{5, 8}, {13, 20}}
```

### Notes

`Dot` contracts the last index of each argument with the first index of the next, so it covers matrix-matrix products, matrix-vector products, and the vector dot product (which yields a scalar) with a single operator. The infix `.` form is equivalent to `Dot[...]` and chains for products of three or more arrays. Exact and symbolic inputs use the elementwise sum of products; machine-precision Real and Complex matrix-matrix products dispatch to a BLAS `gemm` kernel when available.
