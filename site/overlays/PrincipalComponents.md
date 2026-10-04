### Worked examples

```mathematica
In[1]:= PrincipalComponents[{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]  (* collinear points leave the second component at exactly 0 *)
```

```mathematica
In[1]:= p = PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}]  (* rows are observations, columns variables *)
```

```mathematica
In[1]:= {Variance[Map[First, p]], Variance[Map[Last, p]]}  (* variance concentrates in the leading component *)
```

```mathematica
In[1]:= PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}, Method -> "Correlation"]  (* standardise to unit variance first when the units differ *)
```

### Notes

The transform is an orthogonal rotation, so total variance is preserved — it is
redistributed into the leading components, not created or destroyed. Rank-deficient
input therefore puts exactly zero variance in the trailing components, as the collinear
first example shows.

`Method -> "Correlation"` standardises each variable to unit variance first, which is
what you want when the columns have incommensurable units, otherwise the variable with
the largest raw scale dominates for no statistical reason; the default `"Covariance"`
does not. Eigenvalues are not comparable across the two.

Eigenvector signs are canonical — each component is flipped so its largest-magnitude
loading is positive — so the output does not depend on whether the binary was linked
against LAPACK or fell back to the in-house Jacobi solver.
