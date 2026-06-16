---
status: Stable
references:
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — the symmetric and unsymmetric eigenvalue problems."
  - "L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, SIAM, 1997 — eigenvalue algorithms and the QR iteration."
---
### Worked examples

```mathematica
In[1]:= Eigenvalues[{{2, 1}, {0, 3}}]
Out[1]= {3, 2}
```

```mathematica
In[1]:= Eigenvalues[{{2, 0}, {0, 5}}]
Out[1]= {5, 2}
```

```mathematica
In[1]:= Eigenvalues[{{0, 1}, {-1, 0}}]
Out[1]= {-I, I}
```

```mathematica
In[1]:= Eigenvalues[{{a, b}, {c, d}}]
Out[1]= {1/2 (a + d + Sqrt[(-a - d)^2 - 4 (-b c + a d)]), 1/2 (a + d - Sqrt[(-a - d)^2 - 4 (-b c + a d)])}
```

```mathematica
In[1]:= Eigenvalues[{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}}]
Out[1]= {1, -(-1)^(1/3), (-1)^(2/3)}
```

```mathematica
In[1]:= Eigenvalues[{{5, 4, 2}, {4, 5, 2}, {2, 2, 2}}]
Out[1]= {10, 1, 1}
```

### Notes

Eigenvalues of an exact matrix are computed as the roots of the characteristic polynomial `Det[m - lambda I]`, so triangular matrices return their diagonal entries directly and a rotation matrix returns the complex conjugate pair `{-I, I}`. Results are ordered by decreasing absolute value, and repeated eigenvalues appear with their full algebraic multiplicity. Symbolic matrices return closed-form roots (the `2x2` case uses the quadratic formula). Approximate Real / MPFR matrices flow through dedicated numerical kernels (Householder + QR for the symmetric path, Hessenberg + Francis QR for the general path) selectable through the `Method` option.
