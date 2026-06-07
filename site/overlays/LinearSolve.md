---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — solving linear systems by Gaussian elimination."
  - "L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — LU factorisation and triangular solves."
---
### Worked examples

```mathematica
In[1]:= LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]
Out[1]= {-4, 9/2}
```

```mathematica
In[1]:= LinearSolve[{{2, 0}, {0, 4}}, {6, 8}]
Out[1]= {3, 2}
```

```mathematica
In[1]:= LinearSolve[{{1, 1}, {0, 1}}, {3, 1}]
Out[1]= {2, 1}
```

```mathematica
In[1]:= LinearSolve[{{1, 2}, {2, 4}}, {1, 3}]
Out[1]= LinearSolve[{{1, 2}, {2, 4}}, {1, 3}]
```

### Notes

`LinearSolve[m, b]` returns an `x` satisfying `m . x == b`; over exact integer inputs the solution is exact rational, as in the first example. The default algorithm is fraction-free (Bareiss-like) Gauss-Jordan elimination on the augmented matrix `[m | b]`, so no spurious denominators appear during elimination. When the system is inconsistent (last example) `LinearSolve::nosol` is emitted on stderr and the call is returned unevaluated; for under-determined systems a particular solution with free variables set to 0 is returned. The right-hand side may also be a matrix (one column per system), giving a matrix of solutions.
